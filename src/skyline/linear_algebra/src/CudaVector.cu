#include <cassert>

#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaVector.h"

namespace cuda_linalg
{

template<typename T>
CudaVector<T>::CudaVector()
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, int>)
    data_type = CUDA_R_32I;
  else if constexpr (std::is_same_v<T, int64_t>)
    data_type = CUDA_R_64I;
  else if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  //T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  //CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, 0, d_data_ptr, data_type) )
}

template<typename T>
CudaVector<T>::CudaVector(const std::vector<T>& h_vector) 
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, int>)
    data_type = CUDA_R_32I;
  else if constexpr (std::is_same_v<T, int64_t>)
    data_type = CUDA_R_64I;
  else if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  size_t n = h_vector.size();
  d_data_.resize(n);
  thrust::copy(h_vector.begin(), h_vector.end(), d_data_.begin());

  //T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  //CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, data_type) )
}

template<typename T>
CudaVector<T>::CudaVector(size_t n)
{
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, int>)
    data_type = CUDA_R_32I;
  else if constexpr (std::is_same_v<T, int64_t>)
    data_type = CUDA_R_64I;
  else if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  d_data_.resize(n); 
  thrust::fill(d_data_.begin(), d_data_.end(), 0);

  //T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  //CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, data_type) )
}

template<typename T>
void 
CudaVector<T>::resetDescriptor()
{
  CHECK_CUSPARSE(
    cusparseDnVecSetValues(
    descriptor_, 
    thrust::raw_pointer_cast(d_data_.data())) )
}

template<typename T>
void 
CudaVector<T>::resize(size_t n, T val) 
{ 
  cudaDataType_t data_type;
  if constexpr (std::is_same_v<T, int>)
    data_type = CUDA_R_32I;
  else if constexpr (std::is_same_v<T, int64_t>)
    data_type = CUDA_R_64I;
  else if constexpr (std::is_same_v<T, float>)
    data_type = CUDA_R_32F;
  else if constexpr (std::is_same_v<T, double>)
    data_type = CUDA_R_64F;
  else
    assert(0);

  d_data_.clear(); 
  d_data_.resize(n); 
  thrust::fill(d_data_.begin(), d_data_.end(), val);

  //CHECK_CUSPARSE( cusparseDestroyDnVec(descriptor_) )
  //T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  //CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, data_type) )
}

// This is to separate header and .cu
template CudaVector<int>::CudaVector();
template CudaVector<int>::CudaVector(size_t n);
template CudaVector<int>::CudaVector(const std::vector<int>& h_vector);
template void CudaVector<int>::resize(size_t n, int val);
template void CudaVector<int>::resetDescriptor();

template CudaVector<int64_t>::CudaVector();
template CudaVector<int64_t>::CudaVector(size_t n);
template CudaVector<int64_t>::CudaVector(const std::vector<int64_t>& h_vector);
template void CudaVector<int64_t>::resize(size_t n, int64_t val);
template void CudaVector<int64_t>::resetDescriptor();

template CudaVector<float>::CudaVector();
template CudaVector<float>::CudaVector(size_t n);
template CudaVector<float>::CudaVector(const std::vector<float>& h_vector);
template void CudaVector<float>::resize(size_t n, float val);
template void CudaVector<float>::resetDescriptor();

template CudaVector<double>::CudaVector();
template CudaVector<double>::CudaVector(size_t n);
template CudaVector<double>::CudaVector(const std::vector<double>& h_vector);
template void CudaVector<double>::resize(size_t n, double val);
template void CudaVector<double>::resetDescriptor();

}
