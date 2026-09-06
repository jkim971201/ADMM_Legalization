#ifndef DEVICE_UTIL_HPP
#define DEVICE_UTIL_HPP

#include <limits>

#include "GlobalUtil.h"

namespace legalizer
{

__device__ inline bool isPlacedInCorrectRow(
  const int cell_h,
  const int cell_has_ground_at_bottom,
  const int is_row_vdd_up)
{
  // For even-height cell, we have to check for orientation.
  // == has higher priority than bit-wise AND (&).
  if((cell_h & 1) == 0 and cell_has_ground_at_bottom != is_row_vdd_up)
    return false;
  else
    return true;
}

__device__ inline bool isOutOfGrid(
  const int x_size,
  const int y_size,
  const int cell_w,
  const int cell_h,
  const int new_x_grid, 
  const int new_y_grid)
{
  if(new_x_grid + cell_w > x_size || new_y_grid + cell_h > y_size)
    return true;
  else if(new_x_grid < 0 || new_y_grid < 0)
    return true;
  else
    return false;
}

// This must not be called frequently
// (there's quite large overhead)
inline size_t getMaxSharedMemsize()
{
  int device_id;
  cudaGetDevice(&device_id);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device_id);
  return prop.sharedMemPerBlock;
}

inline size_t getMaxSharedMemsizeOptin()
{
  int device_id;
  cudaGetDevice(&device_id);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device_id);
  return prop.sharedMemPerBlockOptin;
}

}

#endif
