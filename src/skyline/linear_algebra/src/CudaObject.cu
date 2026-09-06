#include "cuda_linalg/CudaObject.h"
#include "cuda_linalg/CudaUtil.h"

namespace cuda_linalg
{

// Initialize Static member
CudaHandleSet CudaObject::handle_set = CudaHandleSet();

CudaHandleSet::CudaHandleSet()
{
  CHECK_CUSPARSE( cusparseCreate(&handle_cusparse) );
  CHECK_CUSOLVER( cusolverDnCreate(&handle_cusolver) );
  CHECK_CUBLAS( cublasCreate(&handle_cublas) );
}

}
