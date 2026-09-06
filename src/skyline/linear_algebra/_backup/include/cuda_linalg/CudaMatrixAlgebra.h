#ifndef CUDA_MATRIX_ALGEBRA_H
#define CUDA_MATRIX_ALGEBRA_H

#include "cuda_linalg/CudaMatrix.h"

namespace cuda_linalg
{

/* Compute B = B + alpha * A */
void fmatrixAxpy(const double alpha,
                 const CudaFlattenMatrix<double>& A,
                       CudaFlattenMatrix<double>& B);

/* Compute C = alpha * A + beta * B */
void fmatrixAdd(const double alpha,
                const double beta,
                const CudaFlattenMatrix<double>& A,
                const CudaFlattenMatrix<double>& B,
                      CudaFlattenMatrix<double>& C);

/* Compute C = Tr(A X B) */
double fmatrixInnerProduct(const CudaFlattenMatrix<double>& A,
                           const CudaFlattenMatrix<double>& B);

/* Compute EigenDecomposition A = V * diag(w) * V' */
void computeEVD(const CudaFlattenMatrix<double>& d_A, 
                  CudaVector<int>&    d_eigen_info,
                  CudaVector<double>& d_workspace,
                  CudaVector<double>& d_eigen_values,
                  CudaVector<double>& d_eigen_vectors);

/* Compute Rank-1 Update  A = A + alpha * x * y^T (A is dense matrix) */
void rank1Update(const double alpha,
                 const CudaVector<double>& d_x,
                 const CudaVector<double>& d_y,
                       CudaVector<double>& d_dense_A);

/* Compute Rank-1 Update  A = A + alpha * x * y^T (A is dense matrix) */
void rank1UpdatePtr(const    int  n,  // n : vector size
                    const double  alpha,
                    const double* d_x,
                          CudaVector<double>& d_dense_A);

/* Compute y = alpha * A * x */
void computeMatrixVectorProduct(const double alpha,
                                const CudaFlattenMatrix<double>& A,
                                const CudaVector<double>& x,
                                      CudaVector<double>& y);

/* Compute C = alpha * A * B (A must be symmetric) */
void computeSymMatrixMult(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
  const CudaFlattenMatrix<double>& B,
        CudaFlattenMatrix<double>& C);

/* C = alpha * A * A^T + beta * C */
void symmetricRankKUpdate1(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
        CudaFlattenMatrix<double>& C);

/* C = alpha * A * B^T + alpha * B * A^T + beta * C */
void symmetricRankKUpdate2(
  const double alpha,
  const CudaFlattenMatrix<double>& A,
  const CudaFlattenMatrix<double>& B,
        CudaFlattenMatrix<double>& C);

/* Print Flattened Matrix in Row Major Format */
void printDenseMatrixRowMajor(const CudaFlattenMatrix<double>& d_A,
                              const std::string& title = "");

}

#endif
