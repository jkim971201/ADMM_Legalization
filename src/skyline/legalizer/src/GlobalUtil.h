#ifndef GLOBAL_UTIL_H
#define GLOBAL_UTIL_H

namespace legalizer
{

constexpr float k_infinity = std::numeric_limits<float>::infinity();
constexpr int k_int_max = std::numeric_limits<int>::max();
constexpr int k_int_min = std::numeric_limits<int>::min();

constexpr int k_shape_regular = 0;
constexpr int k_shape_thin_h  = 1;
constexpr int k_shape_thin_v  = 2;
constexpr int k_blockage      = 3;

inline int findMinPowerOfTwoNoLessThan(int n)
{
  int res = 1;
  while(res < n)
    res <<= 1;
  return res;
}

}

#endif
