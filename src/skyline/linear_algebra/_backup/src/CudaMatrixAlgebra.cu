#include <cassert>
#include <vector>

#include <thrust/device_vector.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaMatrixAlgebra.h"
#include "cuda_linalg/CudaVectorAlgebra.h"
#include "cuda_linalg/CudaMatrix.h"

namespace cuda_linalg
{

__global__ void fillLowerPartwithUpperPart(int n, double* d_matrix)
{
	const int row = blockIdx.x;
	const int col = blockIdx.y;

  if(row < n && col < n)
  {
    if(row > col) // lower
    {
      int index1 = row + n * col; // lower index
      int index2 = col + n * row; // upper index (opposite)
      d_matrix[index1] = d_matrix[index2];
    }
  }
}

void fmatrixAxpy(const double alpha,
                 const CudaFlattenMatrix<double>& A,
                       CudaFlattenMatrix<double>& B)
{
  int num_row_A = A.getNumRow();
  int num_col_A = A.getNumCol();

  int num_row_B = B.getNumRow();
  int num_col_B = B.getNumCol();

  assert(num_row_A == num_row_B);
  assert(num_col_A == num_col_B);

  const auto& d_data_A = A.getFlattenVector();
        auto& d_data_B = B.getFlattenVector();

  vectorAxpy(alpha, d_data_A, d_data_B);
}

void fmatrixAdd(const double alpha,
                const double beta,
                const CudaFlattenMatrix<double>& A,
                const CudaFlattenMatrix<double>& B,
                      CudaFlattenMatrix<double>& C)
{
  int num_row_A = A.getNumRow();
  int num_col_A = A.getNumCol();

  int num_row_B = B.getNumRow();
  int num_col_B = B.getNumCol();

  int num_row_C = C.getNumRow();
  int num_col_C = C.getNumCol();

  assert(num_row_A == num_row_B);
  assert(num_row_B == num_row_C);

  assert(num_col_A == num_col_B);
  assert(num_col_B == num_col_C);

  const auto& d_data_A = A.getFlattenVector();
  const auto& d_data_B = B.getFlattenVector();
        auto& d_data_C = C.getFlattenVector();

  vectorAdd(alpha, beta, d_data_A, d_data_B, d_data_C);
}

double fmatrixInnerProduct(const CudaFlattenMatrix<double>& A,
                           const CudaFlattenMatrix<double>& B)
{
  const auto& d_A_data = A.getFlattenVector();
  const auto& d_B_data = B.getFlattenVector();
  return innerProduct(d_A_data, d_B_data);
}

void computeEVD(const CudaFlattenMatrix<double>& d_A,              
                  CudaVector<int>&    d_eigen_info,
                  CudaVector<double>& d_eigen_workspace,
                  CudaVector<double>& d_eigen_values,
                  CudaVector<double>& d_eigen_vectors)
{
  int n = d_A.getNumRow();

  const auto& d_A_data = d_A.getFlattenVector();
  thrust::copy(d_A_data.begin(), d_A_data.end(),
               d_eigen_vectors.begin());

  // syevd returns eigen_vector in the input dense matrix.
  double* d_dense_ptr        = d_eigen_vectors.data();
  double* d_eigen_values_ptr = d_eigen_values.data();

  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; 
  cublasFillMode_t  uplo = CUBLAS_FILL_MODE_UPPER;

  int lwork = 0;            

  CHECK_CUSOLVER(
      cusolverDnDsyevd_bufferSize(
        d_A.getCuSolverHandle(), 
        jobz, 
        uplo, 
        n, 
        d_dense_ptr, 
        n, 
        d_eigen_values_ptr, 
        &lwork));

  d_eigen_workspace.resize(lwork);
  double* d_eigen_workspace_ptr = d_eigen_workspace.data();

  int* d_info_ptr = d_eigen_info.data();

  // Solve A * V = V * diag(eigen_values)
  CHECK_CUSOLVER(
      cusolverDnDsyevd(
        d_A.getCuSolverHandle(), /* handle  : cusolverDnHandle_t                        */
        jobz,                    /* jobz    : EigenSolver option                        */
        uplo,                    /* uplo    : CUBLAS_FILL_MODE_LOWER or *_UPEER         */
        n,                       /* n       : Number of rows (or columns) of A          */
        d_dense_ptr,             /* A       : This will be filled with eigenvectors     */
        n,                       /* lda     : Leading dimension of A, max(lda, n)       */
        d_eigen_values_ptr,      /* W       : EigenValues (array of size n)             */
        d_eigen_workspace_ptr,   /* work    : Workspace   (array of size lwork)         */
        lwork,                   /* Lwork   : Workspace Size (syevd_bufferSize returns) */
        d_info_ptr))             /* devInfo : 0 -> Success, -i -> ith parameter wrong   */
}

void rank1Update(const double alpha,
                 const CudaVector<double>& d_x,
                 const CudaVector<double>& d_y,
                 CudaVector<double>& d_dense_A)
{
  const int n = static_cast<int>(d_x.size());

  /* cublasDger -> A = alpha * x * y^T + A */
  CHECK_CUBLAS(
    cublasDger(d_x.getCuBlasHandle(),
               n,
               n,
               &alpha,           /* alpha    */
               d_x.data(),       /* x_vector */
               1,                /* incx     */
               d_y.data(),       /* y_vector */
               1,                /* incx     */
               d_dense_A.data(), /* matrix A */
               n) )
}

void rank1UpdatePtr(const    int  n,
                    const double  alpha,
                    const double* eigen_vector,
                    CudaVector<double>& d_dense_A)
{
  double* d_dense_A_ptr = d_dense_A.data();

  /* cublasDger -> A = alpha * x * y^T + A */
  CHECK_CUBLAS(
    cublasDger(d_dense_A.getCuBlasHandle(),
               n,
               n,
               &alpha,        /* alpha    */
               eigen_vector,  /* x_vector */
               1,             /* incx     */
               eigen_vector,  /* y_vector */
               1,             /* incx     */
               d_dense_A_ptr, /* matrix A */
               n) )
}

void computeMatrixVectorProduct(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
  const CudaVector<double>& x,
        CudaVector<double>& y)
{
  const double beta = 0.0;
  const int n = static_cast<int>(x.size());

  CHECK_CUBLAS(
    cublasDgemv(
      A.getCuBlasHandle(),
      CUBLAS_OP_N,
      n,
      n,
      &alpha,
      A.getFlattenVector().data(),
      n,
      x.data(),
      1,
      &beta,
      y.data(),
      1) )
}

void computeSymMatrixMult(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
  const CudaFlattenMatrix<double>& B,
        CudaFlattenMatrix<double>& C)
{
  // C = alpha * A * B + beta * C
  const double beta = 0.0;
  const int num_row_A = A.getNumRow();
  const int num_col_A = A.getNumCol();

  const int num_row_B = B.getNumRow();
  const int num_col_B = B.getNumCol();

  const int num_row_C = C.getNumRow();
  const int num_col_C = C.getNumCol();

  const int lda = num_row_A;
  const int ldb = num_row_B;
  const int ldc = num_row_C;

  // A must be symmetric
  assert(num_row_A == num_col_A);

  // C = A * B
  assert(num_row_A == num_row_C);
  assert(num_col_A == num_row_B);
  assert(num_col_B == num_col_C);

  cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
  cublasSideMode_t side = CUBLAS_SIDE_LEFT;

  CHECK_CUBLAS(
    cublasDsymm(
      A.getCuBlasHandle(),
      side,
      uplo,
      num_row_C,
      num_col_C,
      &alpha, A.getFlattenVector().data(), lda,
              B.getFlattenVector().data(), ldb,
      &beta,  C.getFlattenVector().data(), ldc) )
}

void symmetricRankKUpdate1(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
        CudaFlattenMatrix<double>& C)
{
  // C = alpha * A * A^T + beta * C
  const double beta = 0.0;

  const int num_row_A = A.getNumRow();
  const int num_col_A = A.getNumCol();

  const int num_row_C = C.getNumRow();
  const int num_col_C = C.getNumCol();

  assert(num_row_A == num_col_C);
  assert(num_row_C == num_row_C);

  const int lda = num_row_A;
  const int ldc = num_row_C;
  const int n   = num_row_A;
  const int k   = num_col_A;

  cublasFillMode_t  uplo   = CUBLAS_FILL_MODE_UPPER;
  cublasOperation_t transa = CUBLAS_OP_N;

  CHECK_CUBLAS(
    cublasDsyrk(
      A.getCuBlasHandle(),
      uplo,
      transa,
      n,
      k,
      &alpha, A.getFlattenVector().data(), lda,
      &beta,  C.getFlattenVector().data(), ldc) )

  // At this moment, C is not filled fully.
	dim3 grid_size(n,	n, 1);
	dim3 block_size(1, 1);

  //printDenseMatrixRowMajor(C, "C before fill");
  fillLowerPartwithUpperPart<<<grid_size, block_size>>>(n, C.getFlattenVector().data());
  //printDenseMatrixRowMajor(C, "C after  fill");
}

void symmetricRankKUpdate2(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
  const CudaFlattenMatrix<double>& B,
        CudaFlattenMatrix<double>& C)
{
  // C = alpha * A * B^T + alpha * B * A^T + beta * C
  const double beta = 0.0;

  const int num_row_A = A.getNumRow();
  const int num_col_A = A.getNumCol();

  const int num_row_B = B.getNumRow();
  const int num_col_B = B.getNumCol();

  const int num_row_C = C.getNumRow();
  const int num_col_C = C.getNumCol();

  assert(num_row_C == num_col_C);
  assert(num_row_A == num_row_B);
  assert(num_col_A == num_col_B);

  const int lda = num_row_A;
  const int ldb = num_row_B;
  const int ldc = num_row_C;
  const int n   = num_row_A;
  const int k   = num_col_A;

  cublasFillMode_t  uplo   = CUBLAS_FILL_MODE_UPPER;
  cublasOperation_t transa = CUBLAS_OP_N;

  CHECK_CUBLAS(
    cublasDsyr2k(
      A.getCuBlasHandle(),
      uplo,
      transa,
      n,
      k,
      &alpha, A.getFlattenVector().data(), lda,
              B.getFlattenVector().data(), ldb,
      &beta,  C.getFlattenVector().data(), ldc) )

  // At this moment, C is not filled fully.
	dim3 grid_size(n,	n, 1);
	dim3 block_size(1, 1);

  fillLowerPartwithUpperPart<<<grid_size, block_size>>>(n, C.getFlattenVector().data());
}

void printDenseMatrixRowMajor(
  const CudaFlattenMatrix<double>& d_A,
  const std::string& title)
{
  int num_row = d_A.getNumRow();
  int num_col = d_A.getNumCol();

  std::vector<double> h_A_dense(num_row * num_col);

  const auto& d_A_data = d_A.getFlattenVector();

  thrust::copy(d_A_data.begin(), d_A_data.end(), h_A_dense.begin());

  printf("Print Dense Matrix %s\n", title.c_str());
  for(int i = 0; i < num_row; i++) 
  {
    for(int j = 0; j < num_col; j++) 
      printf("%10f ", h_A_dense[i * num_col + j]);
    printf("\n");
  }
  printf("\n");
}

}
