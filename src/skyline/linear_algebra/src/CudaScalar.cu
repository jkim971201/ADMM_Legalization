#include "cuda_linalg/CudaScalar.h"

namespace cuda_linalg
{

template<typename T>
CudaScalar<T>::CudaScalar()
{
  d_data_.resize(1);
  d_ptr_ = d_data_.data();
}

template CudaScalar<int>::CudaScalar();
template CudaScalar<float>::CudaScalar();
template CudaScalar<double>::CudaScalar();

}
