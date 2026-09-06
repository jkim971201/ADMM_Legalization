#ifndef CUDA_VECTOR_H
#define CUDA_VECTOR_H

#include <vector>

#include <cusparse.h>
#include <thrust/device_vector.h>

#include "cuda_linalg/CudaObject.h"

namespace cuda_linalg
{

template<typename T>
class CudaVector : public CudaObject
{
  public:

    CudaVector();
    CudaVector(const std::vector<T>& h_vector);
    CudaVector(size_t n);

    cusparseDnVecDescr_t getDescriptor()       { return descriptor_; }
    cusparseDnVecDescr_t getDescriptor() const { return descriptor_; }

    size_t size() const { return d_data_.size(); }

          T* data()       { return thrust::raw_pointer_cast(d_data_.data()); }
    const T* data() const { return thrust::raw_pointer_cast(d_data_.data()); }

    // thrust::device_vector<T>::iterator -> auto?
          auto begin()       { return d_data_.begin(); }
    const auto begin() const { return d_data_.begin(); }

          auto end()         { return d_data_.end(); }
    const auto end()   const { return d_data_.end(); }

    const thrust::device_vector<T>& getVector() const { return d_data_; }
          thrust::device_vector<T>& getVector()       { return d_data_; }

    void resetDescriptor();

    void swap(CudaVector& vec)
    {
      // swap is allowed between vectors with same size.
      assert(vec.size() == size());
      d_data_.swap(vec.getVector());
      vec.resetDescriptor();
      resetDescriptor();
    }

    void resize(size_t n, T val = 0);

    void clear() { d_data_.clear(); }

    void fillZero() { cudaMemset(data(), 0, sizeof(T) * d_data_.size()); }

    // this = rhs // copy rhs to this
    void operator=(const std::vector<T>& rhs)
    {
      thrust::copy(rhs.begin(), rhs.end(), begin());
    }

  private:

    thrust::device_vector<T> d_data_;
    cusparseDnVecDescr_t descriptor_;
};

}

#endif
