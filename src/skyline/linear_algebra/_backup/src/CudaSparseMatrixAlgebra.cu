#include <vector>

#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaMatrixAlgebra.h"
#include "cuda_linalg/CudaVectorAlgebra.h"
#include "cuda_linalg/CudaSparseMatrixAlgebra.h"

namespace cuda_linalg
{

void sparseMatrixVectorMult(CudaSparseMatrix<double>& d_A,
                      const CudaVector<double>& d_x,
                            CudaVector<double>& d_y)
{
  double alpha = 1.0;
  double beta  = 0.0;

  CHECK_CUSPARSE(
    cusparseSpMV(
      d_A.getCuSparseHandle(),
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha,
      d_A.getSparseDescriptor(),
      d_x.getDescriptor(),
      &beta,
      d_y.getDescriptor(),
      CUDA_R_64F,
      CUSPARSE_SPMV_ALG_DEFAULT,
      d_A.getBufferPtr()) )
}

void makeDenseRowMajorFromSparseCSR(const CudaSparseMatrix<double>& d_A_sparse,
                                                CudaVector<double>& d_A_dense)
{
  int num_nnz = d_A_sparse.getNumNonzero();
  int num_row = d_A_sparse.getNumRow();
  int num_col = d_A_sparse.getNumCol();
  int num_mat_size = num_row * num_col;

  assert(d_A_dense.size() == num_mat_size);

  int leading_dim = num_col; 
  // Leading dimenstion equals to num_col in row-major

  double* d_A_dense_ptr = d_A_dense.data();

  cusparseSpMatDescr_t sparse_descr = d_A_sparse.getSparseDescriptor();

  cusparseDnMatDescr_t dense_descr;
  CHECK_CUSPARSE( 
      cusparseCreateDnMat(
        &dense_descr,           /* dnMatDescr : Dense Matrix Descriptor */
        num_row,  
        num_col, 
        leading_dim,            /* Leading dimension of the matrix */
        (void*)(d_A_dense_ptr), /* Values of the dense matrix      */
        CUDA_R_64F,             /* cudaDataType                    */ 
        CUSPARSE_ORDER_ROW) )

  void* d_buffer = nullptr;
  size_t buffer_size = 0;

  CHECK_CUSPARSE( 
      cusparseSparseToDense_bufferSize(
        d_A_sparse.getCuSparseHandle(), 
        sparse_descr,
        dense_descr,
        CUSPARSE_SPARSETODENSE_ALG_DEFAULT,
        &buffer_size) )

  CHECK_CUDA( cudaMalloc(&d_buffer, buffer_size) )

  CHECK_CUSPARSE( 
      cusparseSparseToDense(
        d_A_sparse.getCuSparseHandle(), 
        sparse_descr,
        dense_descr,
        CUSPARSE_SPARSETODENSE_ALG_DEFAULT,
        d_buffer) )

  CHECK_CUSPARSE( cusparseDestroyDnMat(dense_descr) )

  CHECK_CUDA( cudaFree(d_buffer) )
}

void printCudaSparseMatrixinDense(const CudaSparseMatrix<double>& mat)
{
  int num_row = mat.getNumRow();
  int num_col = mat.getNumCol();
  int mat_size = num_row * num_col;

  CudaVector<double> d_A_dense(mat_size);

  makeDenseRowMajorFromSparseCSR(mat, d_A_dense);

  std::vector<double> h_A_dense(d_A_dense.size());

  thrust::copy(d_A_dense.begin(), d_A_dense.end(),
               h_A_dense.begin());

  printf("Print Dense Matrix\n");
  for(int i = 0; i < num_row; i++) 
  {
    for(int j = 0; j < num_col; j++) 
      printf("%7.2f ", h_A_dense[i * num_col + j]);
    printf("\n");
  }
  printf("\n");
}

void printCudaSparseMatrixinCSR(const CudaSparseMatrix<double>& mat)
{
  std::vector<double> h_nonzero(mat.getNumNonzero());
  std::vector<int>    h_crow(mat.getNumRow() + 1);
  std::vector<int>    h_col(mat.getNumNonzero());

  thrust::copy(mat.getValues().begin(), mat.getValues().end(),
               h_nonzero.begin());

  thrust::copy(mat.getCRows().begin(), mat.getCRows().end(),
               h_crow.begin());

  thrust::copy(mat.getCols().begin(), mat.getCols().end(),
               h_col.begin());

  for(auto& val : h_nonzero)
    printf("Val : %f\n", val);

  for(auto& val : h_crow)
    printf("CRow : %d\n", val);

  for(auto& val : h_col)
    printf("Col : %d\n", val);
}

}
