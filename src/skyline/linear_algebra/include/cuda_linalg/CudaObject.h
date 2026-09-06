#ifndef CUDA_OBJECT_H
#define CUDA_OBJECT_H

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusparse.h>

namespace cuda_linalg
{

struct CudaHandleSet
{
  cusparseHandle_t   handle_cusparse;
  cusolverDnHandle_t handle_cusolver;
  cublasHandle_t     handle_cublas;

  CudaHandleSet();
};

class CudaObject
{
  public:

    CudaObject() {};

    cusparseHandle_t   getCuSparseHandle() const { return handle_set.handle_cusparse; }
    cusolverDnHandle_t getCuSolverHandle() const { return handle_set.handle_cusolver; }
    cublasHandle_t     getCuBlasHandle()   const { return handle_set.handle_cublas;   }

    cusparseHandle_t   getCuSparseHandle() { return handle_set.handle_cusparse; }
    cusolverDnHandle_t getCuSolverHandle() { return handle_set.handle_cusolver; }
    cublasHandle_t     getCuBlasHandle()   { return handle_set.handle_cublas;   }

    static CudaHandleSet handle_set;
};

}

#endif
