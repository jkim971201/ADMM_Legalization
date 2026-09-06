#ifndef CUDA_SPARSE_MATRIX_ALGEBRA_H
#define CUDA_SPARSE_MATRIX_ALGEBRA_H

#include "cuda_linalg/CudaSparseMatrix.h"
#include <thrust/device_vector.h>

namespace cuda_linalg
{

/* Compute y = Ax, A is not const since buffer of A is used in SpMV. */
template<typename T>
void sparseMatrixVectorMult(
  CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y);

/* Compute y = alpha * Ax + beta * y, 
 * A is not const since buffer of A is used in SpMV */
template<typename T>
void sparseMatrixVectorMult(
        T alpha,
        T beta,
        CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y);

/* Compute y = alpha * (A^T)x + beta * y, 
 * A is not const since buffer of A is used in SpMV. 
 * (alpha, beta are host scalar) */
template<typename T>
void sparseMatrixVectorMultTranspose(
        T alpha,
        T beta,
        CudaSparseMatrix<T>& d_A,
  const CudaVector<T>& d_x,
        CudaVector<T>& d_y);

/* Print Sparse Matrix in Dense Matrix Format */
template<typename T>
void printCudaSparseMatrixInDense(const CudaSparseMatrix<T>& mat);

/* Make Dense Matrix (Row Major) from Sprse Matrix (CSR) */
template<typename T>
void makeDenseRowMajorFromSparseCSR(const CudaSparseMatrix<T>& d_A_sparse, CudaVector<T>& d_A_dense);

/* Print Sparse Matrix in CSR Format */
template<typename T>
void printCudaSparseMatrixInCSR(const CudaSparseMatrix<T>& mat);

}

#endif
