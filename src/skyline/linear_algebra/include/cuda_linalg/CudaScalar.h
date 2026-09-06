#ifndef CUDA_SCALAR_H
#define CUDA_SCALAR_H

#include <thrust/device_vector.h>
#include "cuda_linalg/CudaObject.h"

namespace cuda_linalg
{

template<typename T>
class CudaScalar : public CudaObject
{
  public:

    CudaScalar();

          T* data()       { return thrust::raw_pointer_cast(d_ptr_); }
    const T* data() const { return thrust::raw_pointer_cast(d_ptr_); }

          thrust::device_ptr<T> getThrustPtr()       { return d_ptr_; }
    const thrust::device_ptr<T> getThrustPtr() const { return d_ptr_; }

    T getHostValue() const 
    {
      T val;
      cudaMemcpy(&val, data(), sizeof(T), cudaMemcpyDeviceToHost);
      return val;
    }

    void operator=(const T val) { thrust::fill(d_ptr_, d_ptr_ + 1, val); }

  private:

    thrust::device_ptr<T> d_ptr_;
    thrust::device_vector<T> d_data_; 
    // For automatical memory free, we use thrust::device_vector,
    // even for scalar data...
};

}

#endif
