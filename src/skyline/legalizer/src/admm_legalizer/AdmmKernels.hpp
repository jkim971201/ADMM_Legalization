#ifndef ADMM_LEGALIZATION_KERNEL_HPP
#define ADMM_LEGALIZATION_KERNEL_HPP

#include "DeviceUtil.hpp"

#include <thrust/tuple.h>

namespace admm_legalizer
{

struct dual_update_functor
{
  float rho_;
  dual_update_functor(float rho) : rho_(rho) {}

  __host__ __device__ 
  void operator()(thrust::tuple<float&, int> t) 
  {
    float dual = thrust::get<0>(t);
    float usage = float(thrust::get<1>(t));
    float new_dual = dual + rho_ * max(usage - 1, -dual / rho_);
    thrust::get<0>(t) = max(new_dual, 0.0f);
  };
};

struct dual_update_functor_with_step_size
{
  float dual_step_;
  float rho_;
  dual_update_functor_with_step_size(float dual_step, float rho) 
    : dual_step_(dual_step), rho_(rho) {}

  __host__ __device__ 
  void operator()(thrust::tuple<float&, int> t) 
  {
    float lambda = thrust::get<0>(t);
    float usage = float(thrust::get<1>(t));
    float g_tilde = max(usage - 1.0f, -1.0f / rho_);
    float delta = g_tilde + rho_ * 0.5 * g_tilde * g_tilde;
    thrust::get<0>(t) = max(lambda + dual_step_ * delta, 0.0f);
  };
};

__global__ void computeDemandKernel(
  const int    num_cells,
  const int    x_size,
  const int    y_size,
  const int*   lx_grid,
  const int*   ly_grid,
  const int*   w_grid,
  const int*   h_grid,
  const int*   pixel_usage,
        float* cell_demand_plus_area)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_cells)
    return;

  int x_min = lx_grid[thread_id];
  int y_min = ly_grid[thread_id];
  int x_max = x_min + w_grid[thread_id];
  int y_max = y_min + h_grid[thread_id];

  float demand_sum = 0.0f;
  for(int x = x_min; x < x_max; x++)
  {
    if(x < 0 || x >= x_size)
      return;
    for(int y = y_min; y < y_max; y++)
    {
      if(y < 0 || y >= y_size)
        return;

      int pixel_id = x * y_size + y;
      float usage = float(pixel_usage[pixel_id]);
      demand_sum += usage;
    }
  }

  cell_demand_plus_area[thread_id] = demand_sum;
}

__global__ void computeInverseAugmentedDemandKernel(
  const int    num_cells,
  const int    x_size,
  const int    y_size,
  const float  rho,
  const float  expected_avg_disp,
  const int*   lx_grid,
  const int*   ly_grid,
  const int*   w_grid,
  const int*   h_grid,
  const int*   pixel_usage,
        float* cell_inverse_aug_demand)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_cells)
    return;

  int x_min = lx_grid[thread_id];
  int y_min = ly_grid[thread_id];
  int x_max = x_min + w_grid[thread_id];
  int y_max = y_min + h_grid[thread_id];

  float augmented_sum = 0.0f;
  for(int x = x_min; x < x_max; x++)
  {
    if(x < 0 || x >= x_size)
      return;
    for(int y = y_min; y < y_max; y++)
    {
      if(y < 0 || y >= y_size)
        return;
      float usage = float(pixel_usage[x * y_size + y]);
      augmented_sum += 1.0f + usage + rho * 0.5f * usage * usage;
    }
  }

  cell_inverse_aug_demand[thread_id] = expected_avg_disp / augmented_sum;
}

__global__ void assignPartitionKernel(
  const int  num_cells,
  const int  x_hint,
  const int  y_hint,
  const int  x_size,
  const int  y_size,
  const int  offset_x,
  const int  offset_y,
  const int  num_partition_col,
  const int  num_partition_row,
  const int* cell_id_to_w_in_grid,
  const int* cell_id_to_h_in_grid,
  const int* cell_id_to_lx_in_grid,
  const int* cell_id_to_ly_in_grid,
        int* cell_id_to_partition_id)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_cells)
    return;

  const int lx_grid = cell_id_to_lx_in_grid[thread_id] + offset_x;
  const int ly_grid = cell_id_to_ly_in_grid[thread_id] + offset_y;

  const int ux_grid = lx_grid + cell_id_to_w_in_grid[thread_id] - 1;
  const int uy_grid = ly_grid + cell_id_to_h_in_grid[thread_id] - 1;

  const int col1 = (lx_grid / x_hint) % num_partition_col;
  const int col2 = (ux_grid / x_hint) % num_partition_col;

  const int row1 = (ly_grid / y_hint) % num_partition_row;
  const int row2 = (uy_grid / y_hint) % num_partition_row;

  int partition_id = num_partition_col * num_partition_row;
  // default index is an invalid index.
  if(col1 == col2 && row1 == row2)
    partition_id = col1 * num_partition_row + row1; // valid partition index
  // if col1 != col2 or row1 != row2 then
  // we assign invalid partition index to exclude this cell.

  cell_id_to_partition_id[thread_id] = partition_id;
}

__global__ void paintPixelsKernel(
  const int  num_pixels,
  const int  x_hint,
  const int  y_hint,
  const int  x_size,
  const int  y_size,
  const int  offset_x,
  const int  offset_y,
  const int  num_partition_col,
  const int  num_partition_row,
        int* pixel_color)
{
  int pixel_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(pixel_id >= num_pixels)
    return;

  const int pixel_x = pixel_id / y_size + offset_x;
  const int pixel_y = pixel_id % y_size + offset_y;

  const int col = (pixel_x / x_hint) % num_partition_col;
  const int row = (pixel_y / y_hint) % num_partition_row;

  //if(pixel_id / y_size == 0 and (pixel_id % y_size == 1410 or pixel_id % y_size == 1411))
  //  printf("Pixel (%d, %d) offset (%d, %d) row, col (%d, %d) num_col num_row %d %d part_id: %d\n", 
  //    pixel_id / y_size, pixel_id % y_size, offset_x, offset_y, row, col, num_partition_col, num_partition_row, col * num_partition_row + row);

  pixel_color[pixel_id] = col * num_partition_row + row;
}

__device__ inline int computeNumPinVio(
  const bool   will_be_rotated,
  const int    new_x,
  const int    new_y,
  const int    x_size,
  const int    y_size,
  const int    cell_w,
  const int    cell_h,
  const int    site_width,
  const int    row_height,
  const int    num_pin,
  const int    pin_id_offset,
  const float  tech_penalty,
  const int*   pin_id_to_bgn,
  const int*   pin_id_to_end,
  const int*   pin_id_to_layer,
  const int*   pixel_metal_layer)
{
  int pin_short_violation = 0;
  int pin_access_violation = 0;
  for(int i = 0; i < num_pin; i++)
  {
    int pin_id = i + pin_id_offset;
    int pin_bgn = pin_id_to_bgn[pin_id];
    int pin_end = pin_id_to_end[pin_id];
    int pin_layer = pin_id_to_layer[pin_id];

    int pin_bgn_x = pin_bgn / cell_h;
    int pin_bgn_y = pin_bgn % cell_h;

    int pin_end_x = pin_end / cell_h;
    int pin_end_y = pin_end % cell_h;

    if(will_be_rotated == true)
    {
      int orig_bgn_x = pin_bgn_x;
      int orig_bgn_y = pin_bgn_y;

      pin_bgn_x = cell_w - 1 - pin_end_x;
      pin_bgn_y = cell_h - 1 - pin_end_y;

      pin_end_x = cell_w - 1 - orig_bgn_x;
      pin_end_y = cell_h - 1 - orig_bgn_y;
    }

    int x1 = min(pin_bgn_x, pin_end_x) + new_x;
    int y1 = min(pin_bgn_y, pin_end_y) + new_y;
  
    int x2 = max(pin_bgn_x, pin_end_x) + new_x;
    int y2 = max(pin_bgn_y, pin_end_y) + new_y;

    assert(x1 <= x2);
    assert(y1 <= y2);
  
    for(int x = x1; x <= x2; x++)
    {
      for(int y = y1; y <= y2; y++)
      {
        int pixel_layer = pixel_metal_layer[x * y_size + y];
        if(pixel_layer > 0)
        {
          if(pixel_layer == pin_layer)
            pin_short_violation = 1;
          else if(pixel_layer == pin_layer + 1)
            pin_access_violation = 1;
        }
      }
    }
  }

  return pin_short_violation + pin_access_violation;
}

__device__ inline float computeEdgeSpacingCost(
  const bool   will_be_rotated,
  const int    new_x,
  const int    new_y,
  const int    x_size,
  const int    y_size,
  const int    cell_w,
  const int    cell_h,
  const int    cell_ledge_type,
  const int    cell_redge_type,
  const float  edge_spacing_penalty,
  const int*   edge_spacing_rules,
  const int*   pixel_num_ledge_type1,
  const int*   pixel_num_ledge_type2,
  const int*   pixel_num_redge_type1,
  const int*   pixel_num_redge_type2)
{
  int real_ledge_type = will_be_rotated == true ? cell_redge_type : cell_ledge_type;
  int real_redge_type = will_be_rotated == true ? cell_ledge_type : cell_redge_type;

  float cost_on_ledge = 0.0f;
  if(real_ledge_type > 0)
  {
    int spacing_between_redge_type1 = edge_spacing_rules[real_ledge_type + 1];
    int spacing_between_redge_type2 = edge_spacing_rules[real_ledge_type + 2];

    int num_vio_between_redge_type1 = 0;
    int num_vio_between_redge_type2 = 0;

    // Check spacing between right edge type 1
    for(int x = new_x - spacing_between_redge_type1; x < new_x; x++)
      for(int y = new_y; y < new_y + cell_h; y++)
        num_vio_between_redge_type1 += x >= 0 ? pixel_num_redge_type1[x * y_size + y] : 0;

    // Check spacing between right edge type 2
    for(int x = new_x - spacing_between_redge_type2; x < new_x; x++)
      for(int y = new_y; y < new_y + cell_h; y++)
        num_vio_between_redge_type2 += x >= 0 ? pixel_num_redge_type2[x * y_size + y] : 0;

    cost_on_ledge = num_vio_between_redge_type1 + num_vio_between_redge_type2;
  }
  
  float cost_on_redge = 0.0f;
  if(real_redge_type > 0)
  {
    int spacing_between_ledge_type1 = edge_spacing_rules[real_redge_type + 1];
    int spacing_between_ledge_type2 = edge_spacing_rules[real_redge_type + 2];

    int num_vio_between_ledge_type1 = 0;
    int num_vio_between_ledge_type2 = 0;

    // Check spacing between left edge type 1
    for(int x = new_x + cell_w; x < new_x + cell_w + spacing_between_ledge_type1; x++)
      for(int y = new_y; y < new_y + cell_h; y++)
        num_vio_between_ledge_type1 += x < x_size ? pixel_num_ledge_type1[x * y_size + y] : 0;

    // Check spacing between left edge type 2
    for(int x = new_x + cell_w; x < new_x + cell_w + spacing_between_ledge_type2; x++)
      for(int y = new_y; y < new_y + cell_h; y++)
        num_vio_between_ledge_type2 += x < x_size ? pixel_num_ledge_type2[x * y_size + y] : 0;

    cost_on_redge = num_vio_between_ledge_type1 + num_vio_between_ledge_type2;
  }

  return edge_spacing_penalty * (cost_on_ledge + cost_on_redge);
}

__device__ inline float computePenaltySum(
  const bool   is_standard_admm,
  const bool   partition_exclusive,
  const bool   final_refine,
  const int    partition_id,
  const int    cell_id,
  const int    new_x,
  const int    new_y,
  const int    x_size,
  const int    y_size,
  const int    cell_w,
  const int    cell_h,
  const int    cell_group_id,
  const int    cell_has_ground_at_bottom,
  const float  rho,
  const int*   grid_y_to_is_vdd_up,
  const int*   pixel_color,
  const int*   pixel_group_id,
  const int*   pixel_valid, 
  const int*   pixel_usage, 
  const float* pixel_dual)  
{
  if(legalizer::isOutOfGrid(x_size, y_size, cell_w, cell_h, new_x, new_y) == true)
    return legalizer::k_infinity;

  if(legalizer::isPlacedInCorrectRow(cell_h, cell_has_ground_at_bottom, grid_y_to_is_vdd_up[new_y]) == false)
    return legalizer::k_infinity;

  float penalty_sum = 0.0;
  for(int x = new_x; x < new_x + cell_w; x++)
  {
    for(int y = new_y; y < new_y + cell_h; y++)
    {
      int pixel_id = x * y_size + y;
      if(pixel_valid[pixel_id] == 0 || pixel_group_id[pixel_id] != cell_group_id)
        return legalizer::k_infinity;

      int partition_id_this_pixel = pixel_color[pixel_id];
      if(partition_id != partition_id_this_pixel and partition_exclusive == true)
        return legalizer::k_infinity;

      float usage = float(pixel_usage[pixel_id]);
      float dual = pixel_dual[pixel_id];

      if(final_refine == true and usage > 0.0f)
        return legalizer::k_infinity;

      if(is_standard_admm == true)
        penalty_sum += rho * max(usage - 1.0f + dual / rho, 0.0f);
      else
        penalty_sum += rho * max(dual * (usage - 1.0f) + dual / rho, 0.0f);
    }
  }

  if(final_refine == true)
    return 0.0f;
  else
    return penalty_sum;
}

__device__ inline int computeCellOverflow(
  const int  cur_x,
  const int  cur_y,
  const int  cell_w,
  const int  cell_h,
  const int  y_size,
  const int* pixel_usage)
{
  int sum_overflow = 0;
  for(int x = cur_x; x < cur_x + cell_w; x++)
    for(int y = cur_y; y < cur_y + cell_h; y++)
      sum_overflow += max(pixel_usage[x * y_size + y] - 1, 0);
  return sum_overflow;
}

__device__ inline void changeUsage(
  const int  stride,
  const int  thread_id,
  const int  new_x,
  const int  new_y,
  const int  cell_w,
  const int  cell_h,
  const int  y_size,
  const int  delta,
        int* pixel_usage)
{
  int num_pixel_to_change = cell_w * cell_h;
  for(int i = thread_id; i < num_pixel_to_change; i += stride)
  {
    int x = new_x + i / cell_h;
    int y = new_y + i % cell_h;
    atomicAdd(&pixel_usage[x * y_size + y], delta);
  }
}

__device__ inline void changePixelEdgeType(
  const int  stride,
  const int  thread_id,
  const int  new_x,
  const int  new_y,
  const int  cell_w,
  const int  cell_h,
  const int  cell_ledge_type, // real edge_type considering future rotation
  const int  cell_redge_type, // real edge_type considering future rotation
  const int  y_size,
  const int  delta,
        int* pixel_num_ledge_type1,
        int* pixel_num_ledge_type2,
        int* pixel_num_redge_type1,
        int* pixel_num_redge_type2)
{
  int num_pixel_to_change = cell_w * cell_h;
  for(int i = thread_id; i < num_pixel_to_change; i += stride)
  {
    int x = new_x + i / cell_h;
    int y = new_y + i % cell_h;

    if((i / cell_h) == 0)
    {
      if(cell_ledge_type == 1)
        atomicAdd(&pixel_num_ledge_type1[x * y_size + y], delta);
      else if(cell_ledge_type == 2)
        atomicAdd(&pixel_num_ledge_type2[x * y_size + y], delta);
    }
    if((i / cell_h) == (cell_w - 1)) // NOTE: no else-if!, for 1x site width cell
    {
      if(cell_redge_type == 1)
        atomicAdd(&pixel_num_redge_type1[x * y_size + y], delta);
      else if(cell_redge_type == 2)
        atomicAdd(&pixel_num_redge_type2[x * y_size + y], delta);
    }
  }
}

__device__ inline void doCommitOrRipup(
  const bool is_commit, /* true? -> ripup / false? -> commit */
  const bool will_be_rotated,
  const int  stride,
  const int  thread_id,
  const int  cell_lx_in_grid,
  const int  cell_ly_in_grid,
  const int  cell_w,
  const int  cell_h,
  const int  cell_ledge_type,
  const int  cell_redge_type,
  const int  y_size,
        int* pixel_usage,
        int* pixel_num_ledge_type1,
        int* pixel_num_ledge_type2,
        int* pixel_num_redge_type1,
        int* pixel_num_redge_type2)
{
  int delta = is_commit == true ? +1 : -1;

  changeUsage(
    stride, thread_id, cell_lx_in_grid, cell_ly_in_grid, 
    cell_w, cell_h, y_size, delta, 
    pixel_usage);

  __syncthreads();

  int real_ledge_type = will_be_rotated == true ? cell_redge_type : cell_ledge_type;
  int real_redge_type = will_be_rotated == true ? cell_ledge_type : cell_redge_type;

  changePixelEdgeType(
    stride, thread_id, cell_lx_in_grid, cell_ly_in_grid,
    cell_w, cell_h, real_ledge_type, real_redge_type, y_size, delta,
    pixel_num_ledge_type1, 
    pixel_num_ledge_type2,
    pixel_num_redge_type1, 
    pixel_num_redge_type2);

  __syncthreads();
}

__device__ inline float computeDispCost(
  const bool  flag_iccad17,
  const int   new_x,
  const int   new_y,
  const int   core_lx,
  const int   core_ly,
  const int   site_width,
  const int   row_height,
  const int   max_disp_in_row_height,
  const int   original_x_dbu,
  const int   original_y_dbu,
  const int   num_cells,
  const int   num_cells_same_height,
  const float max_disp_coeff)
{
  const float disp_coeff = flag_iccad17 == false ? 1.0f : float(num_cells) / float(num_cells_same_height);
  int new_x_dbu = new_x * site_width + core_lx;
  int new_y_dbu = new_y * row_height + core_ly;
  int disp_dbu = abs(new_x_dbu - original_x_dbu) + abs(new_y_dbu - original_y_dbu);
  float disp_in_site = float(disp_dbu) / float(site_width);
  float disp_threshold_in_site = float(max_disp_in_row_height) * float(row_height) / float(site_width);
  return disp_coeff * float(disp_in_site) + max_disp_coeff * max(float(disp_in_site - disp_threshold_in_site), 0.0f);
}

__device__ inline bool willBeRotated(int is_even_height, int is_row_vdd_up, int cell_has_ground_at_bottom)
{
  if(is_even_height == 1)
    return false;
  bool will_be_rotated = is_row_vdd_up == cell_has_ground_at_bottom ? false : true;
  return will_be_rotated;
}

__device__ inline void getNewCellPosition(
  const int  direction_index,
  const int  x_size,
  const int  y_size,
  const int  cell_w,
  const int  cell_h,
  const int  cur_x,
  const int  cur_y,
  const int* direction,
        int& new_x, 
        int& new_y)
{
  const int epsilon = (cell_h % 2 == 0) ? 2 : 1;
  //const int epsilon = 1;
  new_x = cur_x + direction[2 * direction_index + 0];
  new_y = cur_y + epsilon * direction[2 * direction_index + 1];

  // Clamp
  new_x = max(0, min(x_size - cell_w, new_x));
  new_y = max(0, min(y_size - cell_h, new_y));
}

__device__ inline void getNewCellPositionWithProjection(
  const int  direction_index,
  const int  x_size,
  const int  y_size,
  const int  cell_w,
  const int  cell_h,
  const int  cur_x,
  const int  cur_y,
  const int* direction,
  const int* pixel_valid,
  const int* pixel_region_id,
  const int* region_lx,
  const int* region_ly,
  const int* region_ux,
  const int* region_uy,
        int& new_x, 
        int& new_y)
{
  const int direction_x = direction[2 * direction_index + 0];
  const int direction_y = direction[2 * direction_index + 1];

  const int epsilon = (cell_h % 2 == 0) ? 2 : 1;

  new_x = cur_x + direction_x;
  new_y = cur_y + epsilon * direction_y;

  int region_id = -1;
  bool is_invalid = false;

  auto get_projected_position = [&] (int rect_lx, int rect_ly, int rect_ux, int rect_uy, int& new_x, int& new_y)
  {
      //                        Go Right    / Go Left
    if(abs(direction_x) >= abs(direction_y))
      new_x = direction_x > 0 ? rect_ux + 1 : rect_lx - 1;
    else 
      new_y = direction_y > 0 ? rect_uy + 1 : rect_ly - 1;
      //                        Go Up       / Go Down
  };

  for(int x = new_x; x < new_x + cell_w; x++)
  {
    for(int y = new_y; y < new_y + cell_h; y++)
    {
      // Range check should come first
      if(x >= x_size or y >= y_size or x < 0 or y < 0)
      {
        is_invalid = true;
        break;
      }
      else if(pixel_valid[x * y_size + y] == 0)
      {
        is_invalid = true;
        region_id = pixel_region_id[x * y_size + y];
        break;
      }
    }
    if(is_invalid == true)
      break;
  }

  if(region_id != -1)
  {
    int reg_lx = region_lx[region_id];
    int reg_ly = region_ly[region_id];
    int reg_ux = region_ux[region_id];
    int reg_uy = region_uy[region_id];

    get_projected_position(reg_lx, reg_ly, reg_ux, reg_uy, new_x, new_y);
  }

  // Clamp
  new_x = max(0, min(x_size - cell_w, new_x));
  new_y = max(0, min(y_size - cell_h, new_y));
}

__global__ void doPrimalUpdateKernel(
  const bool   flag_iccad17,
  const bool   is_standard_admm,
  const bool   perturbation_on,
  const int    iter,
  const int    num_cells,
  const int    num_directions,
  const int    num_directions_padded,
  const int    x_size,
  const int    y_size,
  const int    core_lx,
  const int    core_ly,
  const int    site_width,
  const int    row_height,
  const int    max_disp_in_row_height,
  const float  max_disp_coeff,
  const float  tech_penalty,
  const float  edge_spacing_penalty,
  const float  rho,
  const int*   partition_id_to_offset,
  const int*   direction_default,
  const int*   direction_horizontal,
  const int*   direction_vertical,
  const int*   k_height_to_num_cells,
  const int*   pin_id_to_bgn,
  const int*   pin_id_to_end,
  const int*   pin_id_to_layer,
  const int*   macro_id_to_pin_offsets,
  const int*   cell_id_to_macro_id,
  const int*   cell_id_grouped_by_partition,
  const int*   cell_w_in_grid,
  const int*   cell_h_in_grid,
  const int*   cell_original_lx_in_dbu,
  const int*   cell_original_ly_in_dbu,
  const int*   cell_group_id,
  const int*   cell_id_to_has_ground_at_bottom,
  const int*   cell_id_to_ledge_type,
  const int*   cell_id_to_redge_type,
  const int*   edge_spacing_rules,
  const int*   grid_y_to_is_vdd_up,
  const int*   pixel_color,
  const int*   pixel_valid,
  const int*   pixel_group_id,
  const int*   pixel_shape,
  const int*   pixel_power_metal,
  const float* pixel_dual,
        int*   pixel_usage,
        int*   pixel_num_ledge_type1,
        int*   pixel_num_ledge_type2,
        int*   pixel_num_redge_type1,
        int*   pixel_num_redge_type2,
        int*   cell_lx_in_grid,
        int*   cell_ly_in_grid)
{
  // shared_mem size is 2 * num_element
  extern __shared__ float buffer_cost[]; 
  int* buffer_index = (int*)(buffer_cost + num_directions_padded);

  // We should use shared memory for multiple purposes...
  int* buffer_to_store_ovf = buffer_index;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int row_height_in_grid = row_height / site_width;

  const int num_cells_this_partition
    = partition_id_to_offset[partition_id + 1] - partition_id_to_offset[partition_id];

  const int* cells_this_partition
    = cell_id_grouped_by_partition + partition_id_to_offset[partition_id];

  for(int cell_id_in_part = 0; cell_id_in_part < num_cells_this_partition; cell_id_in_part++)
  {
    const int cell_id = cells_this_partition[cell_id_in_part];

    const int cur_x = cell_lx_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];

    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];

    bool is_placed_in_correct_row 
      = legalizer::isPlacedInCorrectRow(cell_h, cell_has_ground_at_bottom, grid_y_to_is_vdd_up[cur_y]);

    // Store overflow to skip zero overflow cell.
    if(thread_id == 0)
    {
      int cell_ovf = computeCellOverflow(cur_x, cur_y, cell_w, cell_h, y_size, pixel_usage);
      buffer_to_store_ovf[0] = cell_ovf;
    }

    __syncthreads();

    if(perturbation_on == false)
    {
      if(buffer_to_store_ovf[0] == 0 and is_placed_in_correct_row == true)
        continue;
    }

    const int* direction = direction_default;

    /*
    const int shape = pixel_shape[cur_x * y_size + cur_y];
    if(shape == k_shape_thin_h)
      direction = direction_horizontal;
    else if(shape == k_shape_thin_v)
      direction = direction_vertical;
      */

    const int is_even_height = cell_h % 2 == 0 ? 1 : 0;

    const int macro_id = cell_id_to_macro_id[cell_id];
    const int pin_id_offset = macro_id_to_pin_offsets[macro_id];
    const int num_pin = macro_id_to_pin_offsets[macro_id + 1] - pin_id_offset;

    const int num_cells_same_height = k_height_to_num_cells[cell_h - 1];

    const int original_x_dbu = cell_original_lx_in_dbu[cell_id];
    const int original_y_dbu = cell_original_ly_in_dbu[cell_id];

    const int group_id_this_cell = cell_group_id[cell_id];

    const int cell_ledge_type = cell_id_to_ledge_type[cell_id];
    const int cell_redge_type = cell_id_to_redge_type[cell_id];

    bool was_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[cur_y], cell_has_ground_at_bottom);

    // Ripup cell
    doCommitOrRipup(false, // false -> ripup
      was_rotated, stride, thread_id, cur_x, cur_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    for(int i = thread_id; i < num_directions; i += stride)
    {
      int new_x, new_y;
      getNewCellPosition(
        i, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, new_x, new_y);

      float disp_cost 
        = computeDispCost(
            flag_iccad17,
            new_x, new_y, core_lx, core_ly, site_width, row_height, 
            max_disp_in_row_height, original_x_dbu, original_y_dbu, 
            num_cells, num_cells_same_height, max_disp_coeff);

      float penalty_cost 
        = computePenaltySum(
            is_standard_admm,  // is_standard_admm
            true,   // partition_exclusive
            false,  // final_refine
            partition_id,
            cell_id, new_x, new_y, x_size, y_size, cell_w, cell_h, 
            group_id_this_cell, cell_has_ground_at_bottom,
            rho, grid_y_to_is_vdd_up, pixel_color, pixel_group_id, pixel_valid, pixel_usage, pixel_dual);

      bool will_be_rotated = (new_y >= y_size or new_y < 0) ? false :
        willBeRotated(is_even_height, grid_y_to_is_vdd_up[new_y], cell_has_ground_at_bottom);

      // If out of die, we don't compute DRV cost
      // (computeNumPinVio requires the cell to be inside the boundary)
      float pin_vio_cost = 0.0f; 
      if(penalty_cost != k_infinity)
      {
        int num_pin_vio_temp = computeNumPinVio(will_be_rotated,
          new_x, new_y, x_size, y_size, cell_w, cell_h, site_width, row_height,
          num_pin, pin_id_offset, tech_penalty,
          pin_id_to_bgn, pin_id_to_end, pin_id_to_layer,
          pixel_power_metal);
        pin_vio_cost = tech_penalty * row_height_in_grid * num_pin_vio_temp;
      }

      float edge_spacing_cost = penalty_cost == k_infinity ? k_infinity
        : computeEdgeSpacingCost(will_be_rotated,
            new_x, new_y, x_size, y_size, cell_w, cell_h,             
            cell_ledge_type, cell_redge_type,
            edge_spacing_penalty,
            edge_spacing_rules,
            pixel_num_ledge_type1,
            pixel_num_ledge_type2,
            pixel_num_redge_type1,
            pixel_num_redge_type2);

      buffer_cost[i] = disp_cost + penalty_cost + pin_vio_cost + edge_spacing_cost;

      //if(cell_id == 370571 and iter == 1)
      //  printf("  CellID: %d CurXY (%d, %d) NewXY (%d, %d) disp_cost: %f penalty_cost: %f\n", cell_id, cur_x, cur_y, new_x, new_y, disp_cost, penalty_cost);
    }

    __syncthreads();

    /* ---- Start Parallel Reduction for Min ---- */
    for(int i = thread_id; i < num_directions_padded; i += stride)
    {
      // Pad infinity to make 2^N array size
      buffer_cost[i] = i < num_directions ? buffer_cost[i] : k_infinity;
      buffer_index[i] = i;
    }

    __syncthreads();

    for(unsigned int s = num_directions_padded / 2; s > 0; s >>= 1)
    {
      if(thread_id < s)
      {
        if(buffer_cost[thread_id] > buffer_cost[thread_id + s])
        {
          buffer_cost[thread_id] = buffer_cost[thread_id + s];
          buffer_index[thread_id] = buffer_index[thread_id + s];
        }
      }
      __syncthreads();
    }
    /* ---- Finish Parallel Reduction for Min ---- */

    float best_cost = buffer_cost[0];
    int best_index = buffer_index[0];

    int new_x, new_y;
    if(best_cost == k_infinity)
    {
      new_x = cur_x;
      new_y = cur_y;
    }
    else
    {
      getNewCellPosition(
        best_index, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, new_x, new_y);
    }

    bool will_be_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[new_y], cell_has_ground_at_bottom);

    // Commit cell
    doCommitOrRipup(true, // true -> commit
      will_be_rotated, stride, thread_id, new_x, new_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    if(thread_id == 0)
    {
      cell_lx_in_grid[cell_id] = new_x;
      cell_ly_in_grid[cell_id] = new_y;
    }

    __syncthreads();
  }
}

__device__ inline void bitonicSortDescending(
  int stride, int thread_id, int n, int* buffer_key, int* buffer_val)
{
  auto swap = [] (int* arr, int i, int j)
  {
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
  };

  auto ascent_swap = [&] (int i, int j)
  {
    if(buffer_key[i] > buffer_key[j])
    {
      swap(buffer_key, i, j);
      swap(buffer_val, i, j);
    }
  };

  auto descent_swap = [&] (int i, int j)
  {
    if(buffer_key[i] < buffer_key[j])
    {
      swap(buffer_key, i, j);
      swap(buffer_val, i, j);
    }
  };

  for(int k = 2; k <= n; k *= 2)
  {
    for(int j = k / 2; j > 0; j /= 2) 
    {
      for(int i = thread_id; i < n; i += stride)
      {
        int ij = i ^ j;
        if(ij > i) 
        {
          if((i & k) == 0) 
            descent_swap(i, ij);
          else 
            ascent_swap(i, ij);
        }
      }
      __syncthreads();
    }
  }
}

__global__ void extractCongestedCellsKernel(
  const int  num_cells,
  const int  y_size,
  const int* pixel_id_to_usage,
  const int* cell_id_to_w,
  const int* cell_id_to_h,
  const int* cell_id_to_lx,
  const int* cell_id_to_ly,
        int* cell_id_to_is_ovf)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_cells)
    return;

  const int cell_id = thread_id;
  const int lx_grid = cell_id_to_lx[cell_id];
  const int ly_grid = cell_id_to_ly[cell_id];

  const int w_grid = cell_id_to_w[cell_id];
  const int h_grid = cell_id_to_h[cell_id];

  const int sum_ovf 
    = computeCellOverflow(lx_grid, ly_grid, w_grid, h_grid, y_size, pixel_id_to_usage);

  cell_id_to_is_ovf[cell_id] = sum_ovf == 0 ? 0 : 1;
}

__global__ void projectionAwareSearchKernel(
  const bool   flag_iccad17,
  const int    num_cells,
  const int    num_ovf_cells,
  const int    num_directions,
  const int    num_directions_padded,
  const int    x_size,
  const int    y_size,
  const int    core_lx,
  const int    core_ly,
  const int    site_width,
  const int    row_height,
  const int    max_disp_in_row_height,
  const float  max_disp_coeff,
  const float  rho,
  const float  coeff_to_neglect_displace,
  const int*   ovf_cells,
  const int*   direction,
  const int*   k_height_to_num_cells,
  const int*   cell_w_in_grid,
  const int*   cell_h_in_grid,
  const int*   cell_original_lx_in_dbu,
  const int*   cell_original_ly_in_dbu,
  const int*   cell_group_id,
  const int*   cell_id_to_has_ground_at_bottom,
  const int*   cell_id_to_ledge_type,
  const int*   cell_id_to_redge_type,
  const int*   grid_y_to_is_vdd_up,
  const int*   region_x_min,
  const int*   region_y_min,
  const int*   region_x_max,
  const int*   region_y_max,
  const int*   pixel_color,
  const int*   pixel_valid,
  const int*   pixel_group_id,
  const int*   pixel_shape,
  const int*   pixel_region_id,
  const float* pixel_dual,
        int*   pixel_usage,
        int*   pixel_num_ledge_type1,
        int*   pixel_num_ledge_type2,
        int*   pixel_num_redge_type1,
        int*   pixel_num_redge_type2,
        int*   cell_lx_in_grid,
        int*   cell_ly_in_grid)
{
  // shared_mem size is 2 * num_element
  extern __shared__ float buffer_cost[]; 
  int* buffer_index = (int*)(buffer_cost + num_directions_padded);

  // We should use shared memory for multiple purposes...
  int* buffer_to_store_ovf = buffer_index;

  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  for(int i = 0; i < num_ovf_cells; i++)
  {
    const int cell_id = ovf_cells[i];
    const int partition_id = -1; // meaningless
    const int cur_x = cell_lx_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];
    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];
    const int is_even_height = cell_h % 2 == 0 ? 1 : 0;
    bool was_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[cur_y], cell_has_ground_at_bottom);

    if(thread_id == 0)
    {
      int cell_ovf = computeCellOverflow(cur_x, cur_y, cell_w, cell_h, y_size, pixel_usage);
      //printf("CellOvf: %d\n", cell_ovf);
      buffer_to_store_ovf[0] = cell_ovf;
    }

    __syncthreads();

    if(buffer_to_store_ovf[0] == 0)
      continue;

    //if(thread_id == 0)
    //  printf("CellID: %d CellW: %d CellH: %d CurXY (%d, %d)\n", cell_id, cell_w, cell_h, cur_x, cur_y);

    const int original_x_dbu = cell_original_lx_in_dbu[cell_id];
    const int original_y_dbu = cell_original_ly_in_dbu[cell_id];

    const int group_id_this_cell = cell_group_id[cell_id];

    const int cell_ledge_type = cell_id_to_ledge_type[cell_id];
    const int cell_redge_type = cell_id_to_redge_type[cell_id];

    const int num_cells_same_height = k_height_to_num_cells[cell_h - 1];

    // Ripup cell
    doCommitOrRipup(false, // false -> ripup
      was_rotated, stride, thread_id, cur_x, cur_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    for(int j = thread_id; j < num_directions; j += stride)
    {
      int new_x, new_y;
      getNewCellPositionWithProjection(
        j, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, pixel_valid, pixel_region_id,
        region_x_min, region_y_min, region_x_max, region_y_max, new_x, new_y);

      float disp_cost 
        = computeDispCost(
            flag_iccad17,
            new_x, new_y, core_lx, core_ly, site_width, row_height, 
            max_disp_in_row_height, original_x_dbu, original_y_dbu, 
            num_cells, num_cells_same_height, max_disp_coeff);

      float penalty_cost 
        = computePenaltySum(
            false,  // is_standard_admm
            false,  // partition_exclusive
            false,  // final_refine
            partition_id,
            cell_id, new_x, new_y, x_size, y_size, cell_w, cell_h, 
            group_id_this_cell, cell_has_ground_at_bottom,
            rho, grid_y_to_is_vdd_up, pixel_color, pixel_group_id, pixel_valid, pixel_usage, pixel_dual);

      buffer_cost[j] = disp_cost / coeff_to_neglect_displace + penalty_cost;

      //printf("  CurXY (%d, %d) NewXY (%d, %d) disp_cost: %f penalty_cost: %f total_cost: %f\n", cur_x, cur_y, new_x, new_y, disp_cost, penalty_cost, buffer_cost[j]);
    }

    __syncthreads();

    /* ---- Start Parallel Reduction for Min ---- */
    for(int j = thread_id; j < num_directions_padded; j += stride)
    {
      // Pad infinity to make 2^N array size
      buffer_cost[j] = j < num_directions ? buffer_cost[j] : k_infinity;
      buffer_index[j] = j;
    }

    __syncthreads();

    for(unsigned int s = num_directions_padded / 2; s > 0; s >>= 1)
    {
      if(thread_id < s)
      {
        if(buffer_cost[thread_id] > buffer_cost[thread_id + s])
        {
          buffer_cost[thread_id] = buffer_cost[thread_id + s];
          buffer_index[thread_id] = buffer_index[thread_id + s];
        }
      }
      __syncthreads();
    }
    /* ---- Finish Parallel Reduction for Min ---- */

    float best_cost = buffer_cost[0];
    int best_index = buffer_index[0];

    int new_x, new_y;
    if(best_cost == k_infinity)
    {
      new_x = cur_x;
      new_y = cur_y;
    }
    else
    {
      getNewCellPositionWithProjection(
        best_index, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, pixel_valid, pixel_region_id,
        region_x_min, region_y_min, region_x_max, region_y_max, new_x, new_y);
    }

    bool will_be_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[new_y], cell_has_ground_at_bottom);

    //if(thread_id == 0)
    //  printf("CellID: %7d CellW: %2d CellH: %1d (%6d, %6d) -> (%6d, %6d)\n", cell_id, cell_w, cell_h, cur_x, cur_y, new_x, new_y);

    // Commit cell
    doCommitOrRipup(true, // true -> commit
      will_be_rotated, stride, thread_id, new_x, new_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    if(thread_id == 0)
    {
      cell_lx_in_grid[cell_id] = new_x;
      cell_ly_in_grid[cell_id] = new_y;
    }

    __syncthreads();
  }
}

__global__ void doFinalRefineKernel(
  const bool   flag_iccad17,
  const int    num_cells,
  const int    num_directions,
  const int    num_directions_padded,
  const int    x_size,
  const int    y_size,
  const int    core_lx,
  const int    core_ly,
  const int    site_width,
  const int    row_height,
  const int    max_disp_in_row_height,
  const float  max_disp_coeff,
  const float  tech_penalty,
  const float  edge_spacing_penalty,
  const float  rho,
  const int*   partition_id_to_offset,
  const int*   direction_default,
  const int*   direction_horizontal,
  const int*   direction_vertical,
  const int*   k_height_to_num_cells,
  const int*   pin_id_to_bgn,
  const int*   pin_id_to_end,
  const int*   pin_id_to_layer,
  const int*   macro_id_to_pin_offsets,
  const int*   cell_id_to_macro_id,
  const int*   cell_id_grouped_by_partition,
  const int*   cell_w_in_grid,
  const int*   cell_h_in_grid,
  const int*   cell_original_lx_in_dbu,
  const int*   cell_original_ly_in_dbu,
  const int*   cell_group_id,
  const int*   cell_id_to_has_ground_at_bottom,
  const int*   cell_id_to_ledge_type,
  const int*   cell_id_to_redge_type,
  const int*   edge_spacing_rules,
  const int*   grid_y_to_is_vdd_up,
  const int*   pixel_color,
  const int*   pixel_valid,
  const int*   pixel_group_id,
  const int*   pixel_shape,
  const int*   pixel_power_metal,
  const float* pixel_dual,
        int*   pixel_usage,
        int*   pixel_num_ledge_type1,
        int*   pixel_num_ledge_type2,
        int*   pixel_num_redge_type1,
        int*   pixel_num_redge_type2,
        int*   cell_lx_in_grid,
        int*   cell_ly_in_grid)
{
  // shared_mem size is 2 * num_element
  extern __shared__ float buffer_cost[]; 
  int* buffer_index = (int*)(buffer_cost + num_directions_padded);

  // We should use shared memory for multiple purposes...
  int* buffer_to_store_height = buffer_index;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int row_height_in_grid = row_height / site_width;

  const int num_cells_this_partition
    = partition_id_to_offset[partition_id + 1] - partition_id_to_offset[partition_id];

  const int* cells_this_partition
    = cell_id_grouped_by_partition + partition_id_to_offset[partition_id];

  for(int cell_id_in_part = 0; cell_id_in_part < num_cells_this_partition; cell_id_in_part++)
  {
    const int cell_id = cells_this_partition[cell_id_in_part];

    const int cur_x = cell_lx_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];

    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];

    bool is_placed_in_correct_row 
      = legalizer::isPlacedInCorrectRow(cell_h, cell_has_ground_at_bottom, grid_y_to_is_vdd_up[cur_y]);

    // Store overflow to skip zero overflow cell.
    if(thread_id == 0)
      buffer_to_store_height[0] = cell_h;

    __syncthreads();

    if(buffer_to_store_height[0] <= 1)
      continue;

    const int* direction = direction_default;

    const int is_even_height = cell_h % 2 == 0 ? 1 : 0;

    const int macro_id = cell_id_to_macro_id[cell_id];
    const int pin_id_offset = macro_id_to_pin_offsets[macro_id];
    const int num_pin = macro_id_to_pin_offsets[macro_id + 1] - pin_id_offset;

    const int num_cells_same_height = k_height_to_num_cells[cell_h - 1];

    const int original_x_dbu = cell_original_lx_in_dbu[cell_id];
    const int original_y_dbu = cell_original_ly_in_dbu[cell_id];

    const int group_id_this_cell = cell_group_id[cell_id];

    const int cell_ledge_type = cell_id_to_ledge_type[cell_id];
    const int cell_redge_type = cell_id_to_redge_type[cell_id];

    bool was_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[cur_y], cell_has_ground_at_bottom);

    // Ripup cell
    doCommitOrRipup(false, // false -> ripup
      was_rotated, stride, thread_id, cur_x, cur_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    for(int i = thread_id; i < num_directions; i += stride)
    {
      int new_x, new_y;
      getNewCellPosition(
        i, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, new_x, new_y);

      float disp_cost 
        = computeDispCost(
            flag_iccad17,
            new_x, new_y, core_lx, core_ly, site_width, row_height, 
            max_disp_in_row_height, original_x_dbu, original_y_dbu, 
            num_cells, num_cells_same_height, max_disp_coeff);

      float penalty_cost 
        = computePenaltySum(
            false, // is_standard_admm
            true,  // partition_exclusive
            true,  // final_refine
            partition_id,
            cell_id, new_x, new_y, x_size, y_size, cell_w, cell_h, 
            group_id_this_cell, cell_has_ground_at_bottom,
            rho, grid_y_to_is_vdd_up, pixel_color, pixel_group_id, pixel_valid, pixel_usage, pixel_dual);

      bool will_be_rotated = (new_y >= y_size or new_y < 0) ? false :
        willBeRotated(is_even_height, grid_y_to_is_vdd_up[new_y], cell_has_ground_at_bottom);

      // If out of die, we don't compute DRV cost
      // (computeNumPinVio requires the cell to be inside the boundary)
      float pin_vio_cost = 0.0f; 
      if(penalty_cost != k_infinity)
      {
        int num_pin_vio_temp = computeNumPinVio(will_be_rotated,
          new_x, new_y, x_size, y_size, cell_w, cell_h, site_width, row_height,
          num_pin, pin_id_offset, tech_penalty,
          pin_id_to_bgn, pin_id_to_end, pin_id_to_layer,
          pixel_power_metal);
        pin_vio_cost = tech_penalty * row_height_in_grid * num_pin_vio_temp;
      }

      float edge_spacing_cost = penalty_cost == k_infinity ? k_infinity
        : computeEdgeSpacingCost(will_be_rotated,
            new_x, new_y, x_size, y_size, cell_w, cell_h,             
            cell_ledge_type, cell_redge_type,
            edge_spacing_penalty,
            edge_spacing_rules,
            pixel_num_ledge_type1,
            pixel_num_ledge_type2,
            pixel_num_redge_type1,
            pixel_num_redge_type2);

      buffer_cost[i] = disp_cost + penalty_cost + pin_vio_cost + edge_spacing_cost;
    }

    __syncthreads();

    /* ---- Start Parallel Reduction for Min ---- */
    for(int i = thread_id; i < num_directions_padded; i += stride)
    {
      // Pad infinity to make 2^N array size
      buffer_cost[i] = i < num_directions ? buffer_cost[i] : k_infinity;
      buffer_index[i] = i;
    }

    __syncthreads();

    for(unsigned int s = num_directions_padded / 2; s > 0; s >>= 1)
    {
      if(thread_id < s)
      {
        if(buffer_cost[thread_id] > buffer_cost[thread_id + s])
        {
          buffer_cost[thread_id] = buffer_cost[thread_id + s];
          buffer_index[thread_id] = buffer_index[thread_id + s];
        }
      }
      __syncthreads();
    }
    /* ---- Finish Parallel Reduction for Min ---- */

    float best_cost = buffer_cost[0];
    int best_index = buffer_index[0];

    int new_x, new_y;
    if(best_cost == k_infinity)
    {
      new_x = cur_x;
      new_y = cur_y;
    }
    else
    {
      getNewCellPosition(
        best_index, x_size, y_size, cell_w, cell_h, cur_x, cur_y, direction, new_x, new_y);
    }

    bool will_be_rotated = willBeRotated(is_even_height, grid_y_to_is_vdd_up[new_y], cell_has_ground_at_bottom);

    // Commit cell
    doCommitOrRipup(true, // true -> commit
      will_be_rotated, stride, thread_id, new_x, new_y, cell_w, cell_h, 
      cell_ledge_type, cell_redge_type,
      y_size, 
      pixel_usage,
      pixel_num_ledge_type1,
      pixel_num_ledge_type2,
      pixel_num_redge_type1,
      pixel_num_redge_type2);

    if(thread_id == 0)
    {
      cell_lx_in_grid[cell_id] = new_x;
      cell_ly_in_grid[cell_id] = new_y;
    }

    __syncthreads();
  }
}

}

#endif
