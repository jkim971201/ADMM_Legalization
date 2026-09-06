#ifndef CUDA_UTIL_H
#define CUDA_UTIL_H

#include <cstdio>
#include <cstdlib>

#define CHECK_CUDA(func)                                                 \
{                                                                        \
  cudaError_t status = (func);                                           \
  if(status != cudaSuccess) {                                            \
    printf("CUDA API failed at line %d (%s) with error: %s (%d)\n",      \
            __LINE__, __FILE__, cudaGetErrorString(status), status);     \
    exit(EXIT_FAILURE);                                                  \
  }                                                                      \
}

#define CHECK_CUSPARSE(func)                                             \
{                                                                        \
  cusparseStatus_t status = (func);                                      \
  if(status != CUSPARSE_STATUS_SUCCESS) {                                \
    printf("cuSPARSE API failed at line %d (%s) with error: %s (%d)\n",  \
            __LINE__, __FILE__, cusparseGetErrorString(status), status); \
    exit(EXIT_FAILURE);                                                  \
  }                                                                      \
}

#define CHECK_CUBLAS(func)                                               \
{                                                                        \
  cublasStatus_t err_ = (func);                                          \
  if (err_ != CUBLAS_STATUS_SUCCESS) {                                   \
    printf("cublas error %d at %s:%d\n", err_, __FILE__, __LINE__);      \
    exit(EXIT_FAILURE);                                                  \
  }                                                                      \
}

#define CHECK_CUSOLVER(func)                                             \
{                                                                        \
  cusolverStatus_t err_ = (func);                                        \
  if (err_ != CUSOLVER_STATUS_SUCCESS) {                                 \
    printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__);    \
    exit(EXIT_FAILURE);                                                  \
  }                                                                      \
}

#endif
