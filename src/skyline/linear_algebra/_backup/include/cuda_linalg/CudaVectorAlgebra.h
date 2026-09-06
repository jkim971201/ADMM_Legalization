#ifndef CUDA_VECTOR_ALGEBRA_H
#define CUDA_VECTOR_ALGEBRA_H

#include "cuda_linalg/CudaVector.h"

namespace cuda_linalg
{

template<typename T>
struct scalarMul
{
  T a;

  scalarMul(T _a) : a(_a) {}

  __host__ __device__ 
  T operator() (const T& x) const
  {
    return a * x;
  }
};

template<typename T>
struct absolute
{
  __host__ __device__
  T operator() (const T& x) const
  {
    if(x >= 0)
      return x;
    else
      return x * -1;
  }
};

template<typename T>
struct axby2
{
  T a;
  T b;

  axby2(T _a, T _b) : a(_a), b(_b) {}

  __host__ __device__ 
  T operator() (const T& x, const T& y) const
  {
    return a * x + b * y * y; 
  }
};

/* Scalar Multiplication  a = k * a */
template<typename T>
void vectorMulScalar(T k, CudaVector<T>& a);

/* Scalar Multiplication  b = k * a */
template<typename T>
void vectorMulScalar(T k, const CudaVector<T>& a, CudaVector<T>& b);

/* Vector Axpy b = b + alpha * a  */
template<typename T>
void vectorAxpy(T alpha, const CudaVector<T>& a, CudaVector<T>& b);

/* Vector Axpy c = alpha * a + beta * b^2 */
template<typename T>
void vectorAxpySquare(
  T alpha, T beta, 
  const CudaVector<T>& a, 
  const CudaVector<T>& b,
        CudaVector<T>& c);

/* Vector Add c = alpha * a + beta * b  */
template<typename T>
void vectorAdd(
  const T alpha, 
  const T beta,
  const CudaVector<T>& a, 
  const CudaVector<T>& b,
        CudaVector<T>& c);

/* Element-wise Multiplication c = a .* b */
template<typename T>
void vectorElementWiseMult(
  const CudaVector<T>& a,
  const CudaVector<T>& b,
        CudaVector<T>& c);

/* Compute a = a / ||a||_2 */
template<typename T>
void normalizeVector(CudaVector<T>& a);

/* Make Absolute vector b = |a| */
template<typename T>
void makeAbsoluteVector(const CudaVector<T>& a, CudaVector<T>& b);

/* Compute ||a||_1 */
template<typename T>
T compute1Norm(const CudaVector<T>& a);

/* Compute ||a||_2 */
template<typename T>
T compute2Norm(const CudaVector<T>& a);

/* Compute ||a - b||_2 */
template<typename T>
T computeDelta2Norm(
  const CudaVector<T>& a,
  const CudaVector<T>& b,
        CudaVector<T>& workspace); // workspace = a - b

/* Compute Sum(a) */
template<typename T>
T computeVectorSum(const CudaVector<T>& a);

/* Compute a^T * b */
template<typename T>
T innerProduct(const CudaVector<T>& a, const CudaVector<T>& b);

/* Find max(a) */
template<typename T>
T computeVectorMax(const CudaVector<T>& a);

template<typename T>
void vectorInit(CudaVector<T>& a, T x);

template<typename T>
void printVector(const CudaVector<T>& d_vector, const std::string& title);

}

#endif
