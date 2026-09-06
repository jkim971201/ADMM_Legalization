#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaVector.h"

namespace cuda_linalg
{

template<typename T>
CudaVector<T>::CudaVector()
{
  T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, 0, d_data_ptr, CUDA_R_64F) )
}

template<typename T>
CudaVector<T>::CudaVector(const std::vector<T>& h_vector) 
{
  size_t n = h_vector.size();
  d_data_.resize(n);
  thrust::copy(h_vector.begin(), h_vector.end(), d_data_.begin());
  T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, CUDA_R_64F) )
}

template<typename T>
CudaVector<T>::CudaVector(size_t n)
{
  d_data_.resize(n); 
  thrust::fill(d_data_.begin(), d_data_.end(), 0);

  T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, CUDA_R_64F) )
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
  d_data_.clear(); 
  d_data_.resize(n); 
  thrust::fill(d_data_.begin(), d_data_.end(), val);

  CHECK_CUSPARSE( cusparseDestroyDnVec(descriptor_) )
  T* d_data_ptr = thrust::raw_pointer_cast(d_data_.data());
  CHECK_CUSPARSE( cusparseCreateDnVec(&descriptor_, n, d_data_ptr, CUDA_R_64F) )
}

// This is to separate header and .cu
template CudaVector<int>::CudaVector();
template CudaVector<int>::CudaVector(size_t n);
template CudaVector<int>::CudaVector(const std::vector<int>& h_vector);
template void CudaVector<int>::resize(size_t n, int val);
template void CudaVector<int>::resetDescriptor();

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
