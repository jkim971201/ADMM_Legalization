#ifndef CUDA_SPARSE_MATRIX_ALGEBRA_H
#define CUDA_SPARSE_MATRIX_ALGEBRA_H

#include "cuda_linalg/CudaSparseMatrix.h"
#include <thrust/device_vector.h>

namespace cuda_linalg
{

/* Compute y = Ax, A is not const since buffer of A is used in SpMV. */
void sparseMatrixVectorMult(CudaSparseMatrix<double>& d_A,
                      const CudaVector<double>& d_x,
                            CudaVector<double>& d_y);

/* Print Sparse Matrix in Dense Matrix Format */
void printCudaSparseMatrixinDense(const CudaSparseMatrix<double>& mat);

/* Make Dense Matrix (Row Major) from Sprse Matrix (CSR) */
void makeDenseRowMajorFromSparseCSR(const CudaSparseMatrix<double>&  d_A_sparse,
                                                 CudaVector<double>& d_A_dense);

/* Print Sparse Matrix in CSR Format */
void printCudaSparseMatrixinCSR(const CudaSparseMatrix<double>& mat);

}

#endif
