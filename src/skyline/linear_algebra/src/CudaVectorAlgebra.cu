#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaVectorAlgebra.h"

#include <cstdio>
#include <vector>
#include <string>
#include <cublas_v2.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>

namespace cuda_linalg
{

template<typename T>
void vectorMulScalar(T k, CudaVector<T>& a)
{
  thrust::transform(a.begin(), a.end(), 
                    a.begin(), 
                    scalarMul<T>(k));
}

template<typename T>
void vectorMulScalar(T k, const CudaVector<T>& a, CudaVector<T>& b)
{
  thrust::transform(
    a.begin(), a.end(), 
    b.begin(), scalarMul<T>(k));
}

template<typename T>
void vectorAxpy(T alpha, const CudaVector<T>& a, CudaVector<T>& b)
{
  const int n = static_cast<int>(a.size());

  const T* d_a_ptr = a.data();
        T* d_b_ptr = b.data();

  if constexpr (std::is_same_v<T, float>)
  {
    cublasSaxpy( 
        a.getCuBlasHandle(), 
        n, 
        &alpha, 
        d_a_ptr, 
        1,  // incr x
        d_b_ptr, 
        1); // incr y
  }
  else if constexpr (std::is_same_v<T, double>)
  {
    cublasDaxpy( 
      a.getCuBlasHandle(), 
      n, 
      &alpha, 
      d_a_ptr, 
      1,  // incr x
      d_b_ptr, 
      1); // incr y
  }
  else
    assert(0);
}

template<typename T>
void vectorAxpy(
  const CudaScalar<T>& d_alpha, 
  const CudaVector<T>& a, 
        CudaVector<T>& b)
{
  const int n = static_cast<int>(a.size());

  const T* d_a_ptr = a.data();
        T* d_b_ptr = b.data();

  if constexpr (std::is_same_v<T, float>)
  {
    cublasSaxpy( 
        a.getCuBlasHandle(), 
        n, 
        d_alpha.data(),
        d_a_ptr, 
        1,  // incr x
        d_b_ptr, 
        1); // incr y
  }
  else if constexpr (std::is_same_v<T, double>)
  {
    cublasDaxpy( 
      a.getCuBlasHandle(), 
      n, 
      d_alpha.data(),
      d_a_ptr, 
      1,  // incr x
      d_b_ptr, 
      1); // incr y
  }
  else
    assert(0);
}

template<typename T>
void vectorAxpySquare(
  T alpha, T beta, 
  const CudaVector<T>& a, 
  const CudaVector<T>& b,
        CudaVector<T>& c)
{
  thrust::transform(
    a.begin(), a.end(),
    b.begin(), 
    c.begin(),
    axby2<T>(alpha, beta) );
}

template<typename T>
void vectorAdd(
  const T alpha, 
  const T beta,
  const CudaVector<T>& a, 
  const CudaVector<T>& b,
        CudaVector<T>& c)
{
  assert(a.size() == b.size());
  assert(b.size() == c.size());

  c.fillZero();
  //vectorInit(c, static_cast<T>(0.0));

  vectorAxpy(alpha, a, c);
  vectorAxpy(beta, b, c);
}

template<typename T>
void vectorAdd(
  const CudaScalar<T>& alpha, 
  const CudaScalar<T>& beta,
  const CudaVector<T>& a, 
  const CudaVector<T>& b,
        CudaVector<T>& c,
        cudaStream_t stream)
{
  assert(a.size() == b.size());
  assert(b.size() == c.size());

  const int n = static_cast<int>(c.size());
  // Zero is okay when using cudaMemset
  cudaMemsetAsync(c.data(), 0, sizeof(T) * n, stream);
  
  vectorAxpy(alpha, a, c);
  vectorAxpy(beta, b, c);
}

template<typename T>
void vectorElementWiseMult(const CudaVector<T>& a,
                           const CudaVector<T>& b,
                                 CudaVector<T>& c)
{
  thrust::transform(a.begin(), a.end(), 
                    b.begin(), 
                    c.begin(), 
                    thrust::multiplies<T>());
}

template<typename T>
void normalizeVector(CudaVector<T>& a)
{
  T norm2 = computeNorm2(a);
  vectorMulScalar(1.0 / norm2, a);
}

template<typename T>
void makeAbsoluteVector(const CudaVector<T>& a,
                              CudaVector<T>& b)
{
  thrust::transform(a.begin(), a.end(),
                    b.begin(),
                    absolute<T>() );
}

template<typename T>
T computeNorm1(const CudaVector<T>& a)
{
  const int n = static_cast<int>(a.size());
  T abs_sum = 0.0;

  if constexpr(std::is_same_v<T, float>)
  {
    CHECK_CUBLAS(
      cublasSasum(
      a.getCuBlasHandle(),
      n,
      a.data(),
      1, // incr x
      &abs_sum) );
  }
  else if constexpr(std::is_same_v<T, double>)
  {
    CHECK_CUBLAS(
      cublasDasum(
      a.getCuBlasHandle(),
      n,
      a.data(),
      1, // incr x
      &abs_sum) );
  }
  else
    assert(0);

  return abs_sum;
}

/* Compute ||a||_2 */
template<typename T>
T computeNorm2(const CudaVector<T>& a)
{
  T sum = 0.0;
  const int n = static_cast<int>(a.size());
  if constexpr(std::is_same_v<T, float>)
    cublasSnrm2(a.getCuBlasHandle(), n, a.data(), 1, &sum);
  else if constexpr(std::is_same_v<T, double>)
    cublasDnrm2(a.getCuBlasHandle(), n, a.data(), 1, &sum);
  else
    assert(0);

  return sum;
}

template<typename T>
void computeNorm2Device(const CudaVector<T>& a, CudaScalar<T>& nrm)
{
  cublasSetPointerMode(a.getCuBlasHandle(), CUBLAS_POINTER_MODE_DEVICE);

  const int n = static_cast<int>(a.size());
  if constexpr(std::is_same_v<T, float>)
    cublasSnrm2(a.getCuBlasHandle(), n, a.data(), 1, nrm.data());
  else if constexpr(std::is_same_v<T, double>)
    cublasDnrm2(a.getCuBlasHandle(), n, a.data(), 1, nrm.data());
  else
    assert(0);

  cublasSetPointerMode(a.getCuBlasHandle(), CUBLAS_POINTER_MODE_HOST);
}

/* Compute ||a - b||_2 */
template<typename T>
T computeDelta2Norm(
  const CudaVector<T>& a,
  const CudaVector<T>& b,
        CudaVector<T>& workspace) // workspace = a - b
{
  vectorAdd(static_cast<T>(1.0), static_cast<T>(-1.0), a, b, workspace);
  T delta_norm2 = computeNorm2(workspace);
  return delta_norm2;
}

template<typename T>
T computeDelta2NormSquare(
  const CudaVector<T>& a,
  const CudaVector<T>& b)
{
  // ||a - b||^2_2 = ||a||^2_2 - 2 * a^T * b + ||b||^2_2
  const T a_norm2 = computeNorm2(a);
  const T b_norm2 = computeNorm2(b);
  const T a_dot_b = innerProduct(a, b);
  return a_norm2 * a_norm2 - T(2.0) * a_dot_b + b_norm2 * b_norm2;
}

/* Compute Sum(a) */
template<typename T>
T computeVectorSum(const CudaVector<T>& a)
{
  return thrust::reduce(a.begin(), a.end());
}

/* Compute a^T * b */
template<typename T>
T innerProduct(const CudaVector<T>& a, const CudaVector<T>& b)
{
  const int n = static_cast<int>(a.size());
  T sum = 0.0;

  if constexpr(std::is_same_v<T, float>)
  {
    CHECK_CUBLAS(
        cublasSdot(
          a.getCuBlasHandle(),
          n,
          a.data(),
          1, /* incr x */
          b.data(),
          1, /* incr y */
          &sum));
  }
  else if constexpr(std::is_same_v<T, double>)
  {
    CHECK_CUBLAS(
        cublasDdot(
          a.getCuBlasHandle(),
          n,
          a.data(),
          1, /* incr x */
          b.data(),
          1, /* incr y */
          &sum));
  }
  else
    assert(0);

  return sum;
}

/* Find max(a) */
template<typename T>
T computeVectorMax(const CudaVector<T>& a)
{
  return *thrust::max_element(thrust::device, a.begin(), a.end());
}

/* Find min(a) */
template<typename T>
T computeVectorMin(const CudaVector<T>& a)
{
  return *thrust::min_element(thrust::device, a.begin(), a.end());
}

template<typename T>
void vectorInit(CudaVector<T>& a, T x)
{
  thrust::fill(a.begin(), a.end(), x);
}

template<typename T>
void printVector(const CudaVector<T>& d_vector, 
                 const std::string& title)
{
  size_t len = d_vector.size();
  std::vector<T> h_vector(len);
  thrust::copy(d_vector.begin(), d_vector.end(),
               h_vector.begin());

  printf("Vector %s\n", title.c_str());
  if constexpr (std::is_same_v<T, int>)
  {
    for(const auto& val : h_vector)
      printf("%d ", val);
  }
  else
  {
    for(const auto& val : h_vector)
      printf("%f ", val);
  }
  printf("\n");
}

// This is to separate header and .cu
template void vectorMulScalar(float k, CudaVector<float>& a);
template void vectorMulScalar(double k, CudaVector<double>& a);
template void vectorMulScalar(float k, const CudaVector<float>& a, CudaVector<float>& b);
template void vectorMulScalar(double k, const CudaVector<double>& a, CudaVector<double>& b);

template float computeDelta2Norm(
  const CudaVector<float>& a, 
  const CudaVector<float>& b,
        CudaVector<float>& workspace);

template double computeDelta2Norm(
  const CudaVector<double>& a, 
  const CudaVector<double>& b,
        CudaVector<double>& workspace);

template float computeDelta2NormSquare(
  const CudaVector<float>& a, 
  const CudaVector<float>& b);

template double computeDelta2NormSquare(
  const CudaVector<double>& a, 
  const CudaVector<double>& b);

template float  computeNorm1(const CudaVector<float>& a);
template double computeNorm1(const CudaVector<double>& a);

template float  computeNorm2(const CudaVector<float>& a);
template double computeNorm2(const CudaVector<double>& a);

template void computeNorm2Device(const CudaVector<float>& a, CudaScalar<float>& nrm);
template void computeNorm2Device(const CudaVector<double>& a, CudaScalar<double>& nrm);

template void vectorAxpy(
  float alpha, 
  const CudaVector<float>& a, 
        CudaVector<float>& b);

template void vectorAxpy(
  double alpha, 
  const CudaVector<double>& a, 
        CudaVector<double>& b);

template void vectorAxpy(
  const CudaScalar<float>& alpha, 
  const CudaVector<float>& a, 
        CudaVector<float>& b);

template void vectorAxpy(
  const CudaScalar<double>& alpha, 
  const CudaVector<double>& a, 
        CudaVector<double>& b);

template void vectorAxpySquare(
  float alpha, float beta, 
  const CudaVector<float>& a, 
  const CudaVector<float>& b,
        CudaVector<float>& c);

template void vectorAxpySquare(
  double alpha, double beta, 
  const CudaVector<double>& a, 
  const CudaVector<double>& b,
        CudaVector<double>& c);

template void vectorAdd(
  const float alpha, 
  const float beta,
  const CudaVector<float>& a, 
  const CudaVector<float>& b,
        CudaVector<float>& c);

template void vectorAdd(
  const double alpha, 
  const double beta,
  const CudaVector<double>& a, 
  const CudaVector<double>& b,
        CudaVector<double>& c);

template void vectorAdd(
  const CudaScalar<float>& alpha, 
  const CudaScalar<float>& beta,
  const CudaVector<float>& a, 
  const CudaVector<float>& b,
        CudaVector<float>& c,
        cudaStream_t cublas_stream);

template void vectorAdd(
  const CudaScalar<double>& alpha, 
  const CudaScalar<double>& beta,
  const CudaVector<double>& a, 
  const CudaVector<double>& b,
        CudaVector<double>& c,
        cudaStream_t cublas_stream);

template void vectorElementWiseMult(
  const CudaVector<float>& a,
  const CudaVector<float>& b,
        CudaVector<float>& c);

template void vectorElementWiseMult(
  const CudaVector<double>& a,
  const CudaVector<double>& b,
        CudaVector<double>& c);

template float  innerProduct(const CudaVector<float>& a, const CudaVector<float>& b);
template double innerProduct(const CudaVector<double>& a, const CudaVector<double>& b);

template int    computeVectorSum(const CudaVector<int>& a);
template float  computeVectorSum(const CudaVector<float>& a);
template double computeVectorSum(const CudaVector<double>& a);

template int    computeVectorMax(const CudaVector<int>& a);
template float  computeVectorMax(const CudaVector<float>& a);
template double computeVectorMax(const CudaVector<double>& a);

template int    computeVectorMin(const CudaVector<int>& a);
template float  computeVectorMin(const CudaVector<float>& a);
template double computeVectorMin(const CudaVector<double>& a);

template void printVector(const CudaVector<int>& d_vector, const std::string& title);
template void printVector(const CudaVector<float>& d_vector, const std::string& title);
template void printVector(const CudaVector<double>& d_vector, const std::string& title);

}
