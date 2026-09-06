#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaVectorAlgebra.h"

#include <cstdio>
#include <vector>
#include <string>
#include <string_view>
#include <cublas_v2.h>
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

  vectorInit(c, static_cast<T>(0.0));

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
  T norm2 = compute2Norm(a);
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
T compute1Norm(const CudaVector<T>& a)
{
  T sum_abs
    = thrust::transform_reduce(a.begin(), a.end(), 
                               absolute<T>(), 
                               static_cast<T>(0.0),
                               thrust::plus<T>());
  return sum_abs;
}

/* Compute ||a||_2 */
template<typename T>
T compute2Norm(const CudaVector<T>& a)
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

/* Compute ||a - b||_2 */
template<typename T>
T computeDelta2Norm(
  const CudaVector<T>& a,
  const CudaVector<T>& b,
        CudaVector<T>& workspace) // workspace = a - b
{
  vectorAdd(static_cast<T>(1.0), static_cast<T>(-1.0), a, b, workspace);
  T delta_norm2 = compute2Norm(workspace);
  return delta_norm2;
}

/* Compute Sum(a) */
template<typename T>
T computeVectorSum(const CudaVector<T>& a)
{
  const int n = static_cast<int>(a.size());
  T sum = 0.0;

  if constexpr(std::is_same_v<T, float>)
  {
    CHECK_CUBLAS(
      cublasSasum(
      a.getCuBlasHandle(),
      n,
      a.data(),
      1, // incr x
      &sum) );
  }
  else if constexpr(std::is_same_v<T, double>)
  {
    CHECK_CUBLAS(
      cublasDasum(
      a.getCuBlasHandle(),
      n,
      a.data(),
      1, // incr x
      &sum) );
  }
  else
    assert(0);

  return sum;
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
  for(const auto& val : h_vector)
    printf("%f ", val);
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

template float  compute2Norm(const CudaVector<float>& a);
template double compute2Norm(const CudaVector<double>& a);

template float  compute1Norm(const CudaVector<float>& a);
template double compute1Norm(const CudaVector<double>& a);

template void vectorAxpy(
  float alpha, 
  const CudaVector<float>& a, 
        CudaVector<float>& b);

template void vectorAxpy(
  double alpha, 
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

template float  innerProduct(const CudaVector<float>& a, const CudaVector<float>& b);
template double innerProduct(const CudaVector<double>& a, const CudaVector<double>& b);

template float  computeVectorSum(const CudaVector<float>& a);
template double computeVectorSum(const CudaVector<double>& a);

template float  computeVectorMax(const CudaVector<float>& a);
template double computeVectorMax(const CudaVector<double>& a);

}
