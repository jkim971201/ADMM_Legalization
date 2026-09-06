#include <vector>

#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaMatrixAlgebra.h"
#include "cuda_linalg/CudaVectorAlgebra.h"
#include "cuda_linalg/CudaSparseMatrixAlgebra.h"

namespace cuda_linalg
{

template<typename T>
void sparseMatrixVectorMult(
  CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y)
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  T alpha = 1.0;
  T beta  = 0.0;

  CHECK_CUSPARSE(
    cusparseSpMV(
      d_A.getCuSparseHandle(),
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha,
      d_A.getSparseDescriptor(),
      d_x.getDescriptor(),
      &beta,
      d_y.getDescriptor(),
      data_type,
      CUSPARSE_SPMV_ALG_DEFAULT,
      d_A.getBufferPtr()) )
}

template<typename T>
void sparseMatrixVectorMult(
        T alpha,
        T beta,
        CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y)
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  CHECK_CUSPARSE(
    cusparseSpMV(
      d_A.getCuSparseHandle(),
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha,
      d_A.getSparseDescriptor(),
      d_x.getDescriptor(),
      &beta,
      d_y.getDescriptor(),
      data_type,
      CUSPARSE_SPMV_ALG_DEFAULT,
      d_A.getBufferPtr()) )
}

template<typename T>
void sparseMatrixVectorMultTranspose(
        T alpha,
        T beta,
        CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y)
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  CHECK_CUSPARSE(
    cusparseSpMV(
      d_A.getCuSparseHandle(),
      CUSPARSE_OPERATION_TRANSPOSE,
      &alpha,
      d_A.getSparseDescriptor(),
      d_x.getDescriptor(),
      &beta,
      d_y.getDescriptor(),
      data_type,
      CUSPARSE_SPMV_ALG_DEFAULT,
      d_A.getBufferPtr()) )
}

template<typename T>
void makeDenseRowMajorFromSparseCSR(
  const CudaSparseMatrix<T>& d_A_sparse,
              CudaVector<T>& d_A_dense)
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  int num_nnz = d_A_sparse.getNumNonzero();
  int num_row = d_A_sparse.getNumRow();
  int num_col = d_A_sparse.getNumCol();
  int num_mat_size = num_row * num_col;

  assert(d_A_dense.size() == num_mat_size);

  int leading_dim = num_col; 
  // Leading dimenstion equals to num_col in row-major

  T* d_A_dense_ptr = d_A_dense.data();

  cusparseSpMatDescr_t sparse_descr = d_A_sparse.getSparseDescriptor();

  cusparseDnMatDescr_t dense_descr;
  CHECK_CUSPARSE( 
      cusparseCreateDnMat(
        &dense_descr,           /* dnMatDescr : Dense Matrix Descriptor */
        num_row,  
        num_col, 
        leading_dim,            /* Leading dimension of the matrix */
        (void*)(d_A_dense_ptr), /* Values of the dense matrix      */
        data_type,              /* cudaDataType                    */ 
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

template<typename T>
void printCudaSparseMatrixInDense(const CudaSparseMatrix<T>& mat)
{
  int num_row = mat.getNumRow();
  int num_col = mat.getNumCol();
  int mat_size = num_row * num_col;

  CudaVector<T> d_A_dense(mat_size);

  makeDenseRowMajorFromSparseCSR(mat, d_A_dense);

  std::vector<T> h_A_dense(d_A_dense.size());

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

template<typename T>
void printCudaSparseMatrixInCSR(const CudaSparseMatrix<T>& mat)
{
  std::vector<T>   h_nonzero(mat.getNumNonzero());
  std::vector<int> h_crow(mat.getNumRow() + 1);
  std::vector<int> h_col(mat.getNumNonzero());

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

// This is to separate header and .cu
template void sparseMatrixVectorMult(
  CudaSparseMatrix<float>& d_A,
  const CudaVector<float>& d_x,
        CudaVector<float>& d_y);

template void sparseMatrixVectorMult(
  CudaSparseMatrix<double>& d_A,
  const CudaVector<double>& d_x,
        CudaVector<double>& d_y);

template void sparseMatrixVectorMult(
  float alpha,
  float beta,
  CudaSparseMatrix<float>& d_A,
  const CudaVector<float>& d_x,
        CudaVector<float>& d_y);

template void sparseMatrixVectorMult(
  double alpha,
  double beta,
  CudaSparseMatrix<double>& d_A,
  const CudaVector<double>& d_x,
        CudaVector<double>& d_y);

template void sparseMatrixVectorMultTranspose(
  float alpha,
  float beta,
  CudaSparseMatrix<float>& d_A,
  const CudaVector<float>& d_x,
        CudaVector<float>& d_y);

template void sparseMatrixVectorMultTranspose(
  double alpha,
  double beta,
  CudaSparseMatrix<double>& d_A,
  const CudaVector<double>& d_x,
        CudaVector<double>& d_y);

template void makeDenseRowMajorFromSparseCSR(
  const CudaSparseMatrix<float>& d_A_sparse,
              CudaVector<float>& d_A_dense);

template void makeDenseRowMajorFromSparseCSR(
  const CudaSparseMatrix<double>& d_A_sparse,
              CudaVector<double>& d_A_dense);

template void printCudaSparseMatrixInDense(const CudaSparseMatrix<float>& mat);
template void printCudaSparseMatrixInDense(const CudaSparseMatrix<double>& mat);

template void printCudaSparseMatrixInCSR(const CudaSparseMatrix<float>& mat);
template void printCudaSparseMatrixInCSR(const CudaSparseMatrix<double>& mat);

}
