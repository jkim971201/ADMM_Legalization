#include <thrust/tuple.h>
#include "DeviceUtil.hpp"

namespace legalizer
{

struct dbu_to_grid_functor
{
  int offset_;
  int unit_;

  dbu_to_grid_functor(int offset, int unit) 
    : offset_(offset), unit_(unit) {}

  __host__ __device__
  int operator()(int dbu) const
  {
    return (dbu - offset_) / unit_;
  }
};

struct grid_to_dbu_functor
{
  int offset_;
  int unit_;

  grid_to_dbu_functor(int offset, int unit)
    : offset_(offset), unit_(unit) {}

  __host__ __device__
  int operator()(int grid) const
  {
    return grid * unit_ + offset_;
  }
};

struct overflow_functor
{
  __host__ __device__ 
  int operator()(int usage) const 
  { 
    int ovf = max(usage - 1, 0);
    return ovf;
  }
};

struct displace_functor
{
  int core_lx_;
  int core_ly_;
  int site_width_;
  int row_height_;
  displace_functor(int core_lx, int core_ly, int site_width, int row_height)
    : core_lx_(core_lx), core_ly_(core_ly), site_width_(site_width), row_height_(row_height) {}

  __host__ __device__  // tuple: {grid_index, dbu}
  float operator()(thrust::tuple<int, int, int, int> t) 
  {
    int grid_x = thrust::get<0>(t);
    int grid_y = thrust::get<1>(t);
    int dbu_x = grid_x * site_width_ + core_lx_;
    int dbu_y = grid_y * row_height_ + core_ly_;
    int original_dbu_x = thrust::get<2>(t);
    int original_dbu_y = thrust::get<3>(t);
    return float(abs(dbu_x - original_dbu_x) + abs(dbu_y - original_dbu_y));
  }
};

__global__ void computeUsageGridKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_w_grid,
  const int* cell_h_grid,
  const int* cell_lx_grid,
  const int* cell_ly_grid,
        int* pixel_usage)
{
  const int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id < num_cells)
  {
    int cell_id = thread_id;

    int w_grid = cell_w_grid[cell_id];
    int h_grid = cell_h_grid[cell_id];

    int lx_grid = cell_lx_grid[cell_id];
    int ly_grid = cell_ly_grid[cell_id];

    int ux_grid = lx_grid + w_grid;
    int uy_grid = ly_grid + h_grid;

    int min_x = max(0, lx_grid);
    int min_y = max(0, ly_grid);

    int max_x = min(x_size, ux_grid);
    int max_y = min(y_size, uy_grid);

    for(int x = min_x; x < max_x; x++)
      for(int y = min_y; y < max_y; y++)
        atomicAdd(&(pixel_usage[x * y_size + y]), 1);
  }
}

__global__ void computePixelEdgeMapKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_id_to_ledge_type,
  const int* cell_id_to_redge_type,
  const int* cell_id_to_has_ground_at_bottom,
  const int* grid_y_to_is_vdd_up,
        int* pixel_num_ledge_type1,
        int* pixel_num_ledge_type2,
        int* pixel_num_redge_type1,
        int* pixel_num_redge_type2)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int x_min = cell_lx_in_grid[cell_id];
    const int y_min = cell_ly_in_grid[cell_id];
    const int x_max = x_min + cell_w;
    const int y_max = y_min + cell_h;

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];

    const int is_row_vdd_up = grid_y_to_is_vdd_up[y_min];

    bool will_be_rotated 
      = is_row_vdd_up == cell_has_ground_at_bottom ? false : true;

    const int cell_ledge_type = cell_id_to_ledge_type[cell_id];
    const int cell_redge_type = cell_id_to_redge_type[cell_id];

    int real_ledge_type = will_be_rotated == true ? cell_redge_type : cell_ledge_type;
    int real_redge_type = will_be_rotated == true ? cell_ledge_type : cell_redge_type;

    for(int x = x_min; x < x_max; x++)
    {
      for(int y = y_min; y < y_max; y++)
      {
        if(x == x_min)
        {
          if(real_ledge_type == 1)
            atomicAdd(&pixel_num_ledge_type1[x * y_size + y], 1);
          else if(real_ledge_type == 2)
            atomicAdd(&pixel_num_ledge_type2[x * y_size + y], 1);
        }
        if(x == x_max - 1) // NOTE: no else-if!, for 1X site width cell
        {
          if(real_redge_type == 1)
            atomicAdd(&pixel_num_redge_type1[x * y_size + y], 1);
          else if(real_redge_type == 2)
            atomicAdd(&pixel_num_redge_type2[x * y_size + y], 1);
        }
      }
    }
  }
}

__device__ inline int getGridStart(int unit, int offset, int dbu)
{
  return (dbu - offset) / unit;
}

__device__ inline int getGridEnd(int unit, int offset, int dbu)
{
  int quo = (dbu - offset) / unit;
  int rem = (dbu - offset) % unit;
  return rem == 0 ? quo - 1 : quo;
}

__global__ void computeUsageDbuKernel(
  const int  num_cells,
  const int  core_lx,
  const int  core_ly,
  const int  site_width,
  const int  row_height,
  const int  x_size,
  const int  y_size,
  const int* cell_w_grid,
  const int* cell_h_grid,
  const int* cell_lx_dbu,
  const int* cell_ly_dbu,
        int* pixel_usage)
{
  const int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id < num_cells)
  {
    int cell_id = thread_id;
    int w_grid = cell_w_grid[cell_id];
    int h_grid = cell_h_grid[cell_id];

    int w_dbu = w_grid * site_width;
    int h_dbu = h_grid * row_height;

    int lx_dbu = cell_lx_dbu[cell_id];
    int ly_dbu = cell_ly_dbu[cell_id];

    int ux_dbu = lx_dbu + w_dbu;
    int uy_dbu = ly_dbu + h_dbu;

    int min_x = getGridStart(site_width, core_lx, lx_dbu);
    int min_y = getGridStart(row_height, core_ly, ly_dbu);

    int max_x = getGridEnd(site_width, core_lx, ux_dbu);
    int max_y = getGridEnd(row_height, core_ly, uy_dbu);

    // Clipping
    min_x = max(min_x, 0);
    min_y = max(min_y, 0);

    max_x = min(max_x, x_size - 1);
    max_y = min(max_y, y_size - 1);

    for(int x = min_x; x <= max_x; x++)
    {
      for(int y= min_y; y <= max_y; y++)
        atomicAdd(&(pixel_usage[x * y_size + y]), 1);
    }
  }
}

__global__ void markPinLayerToGridKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_h_in_grid,
  const int* cell_id_to_cell_macro_id,
  const int* cell_macro_id_to_pin_id_offset,
  const int* pin_id_to_bgn,
  const int* pin_id_to_end,
  const int* pin_id_to_layer,
        int* pixel_pin_layer)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    int lx_grid = cell_lx_in_grid[cell_id];
    int ly_grid = cell_ly_in_grid[cell_id];

    int h_grid = cell_h_in_grid[cell_id];

    int cell_macro_id = cell_id_to_cell_macro_id[cell_id];

    int pin_id_offset = cell_macro_id_to_pin_id_offset[cell_macro_id];
    int num_pin = cell_macro_id_to_pin_id_offset[cell_macro_id + 1] - pin_id_offset;

    for(int i = 0; i < num_pin; i++)
    {
      int pin_id = i + pin_id_offset;
      int pin_bgn = pin_id_to_bgn[pin_id];
      int pin_end = pin_id_to_end[pin_id];
      int pin_layer = pin_id_to_layer[pin_id];

      int x1 = pin_bgn / h_grid + lx_grid;
      int y1 = pin_bgn % h_grid + ly_grid;
  
      int x2 = pin_end / h_grid + lx_grid;
      int y2 = pin_end % h_grid + ly_grid;
  
      for(int x = x1; x <= x2; x++)
      {
        for(int y = y1; y <= y2; y++)
          pixel_pin_layer[x * y_size + y] = pin_layer;
      }
    }
  }
}

__global__ void countPinVioKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_id_to_cell_macro_id,
  const int* cell_macro_id_to_pin_id_offset,
  const int* cell_id_to_has_ground_at_bottom,
  const int* pin_id_to_bgn,
  const int* pin_id_to_end,
  const int* pin_id_to_layer,
  const int* pixel_id_to_layer,
  const int* grid_y_to_is_vdd_up,
        int* cell_id_to_num_short_vio,
        int* cell_id_to_num_access_vio)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int cur_x = cell_lx_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];

    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];

    const int is_row_vdd_up = grid_y_to_is_vdd_up[cur_y];

    const bool will_be_rotated = is_row_vdd_up == cell_has_ground_at_bottom ? false : true;

    const int cell_macro_id = cell_id_to_cell_macro_id[cell_id];

    const int pin_id_offset = cell_macro_id_to_pin_id_offset[cell_macro_id];
    const int num_pin = cell_macro_id_to_pin_id_offset[cell_macro_id + 1] - pin_id_offset;

    bool count_short = false;
    bool count_access = false;
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

      int x1 = min(pin_bgn_x, pin_end_x) + cur_x;
      int y1 = min(pin_bgn_y, pin_end_y) + cur_y;
    
      int x2 = max(pin_bgn_x, pin_end_x) + cur_x;
      int y2 = max(pin_bgn_y, pin_end_y) + cur_y;
  
      assert(x1 <= x2);
      assert(y1 <= y2);

      for(int x = x1; x <= x2; x++)
      {
        for(int y = y1; y <= y2; y++)
        {
          int pixel_layer = pixel_id_to_layer[x * y_size + y];
          if(pixel_layer > 0)
          {
            if(pixel_layer == pin_layer and count_short == false)
            {
              cell_id_to_num_short_vio[cell_id] += 1;
              count_short = true;
            }
            else if(pixel_layer == pin_layer + 1 and count_access == false)
            {
              cell_id_to_num_access_vio[cell_id] += 1;
              count_access = true;
            }
          }
        }
      }
    }
  }
}

__global__ void countEdgeSpacingVioKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_id_to_ledge_type,
  const int* cell_id_to_redge_type,
  const int* cell_id_to_has_ground_at_bottom,
  const int* grid_y_to_is_vdd_up,
  const int* pixel_num_ledge_type1,
  const int* pixel_num_ledge_type2,
  const int* pixel_num_redge_type1,
  const int* pixel_num_redge_type2,
  const int* edge_spacing_rules,
        int* cell_id_to_num_edge_spacing_vio)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int cur_x = cell_lx_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];
    const int cell_w = cell_w_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];

    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];

    const int is_row_vdd_up = grid_y_to_is_vdd_up[cur_y];

    bool will_be_rotated = is_row_vdd_up == cell_has_ground_at_bottom ? false : true;

    const int cell_ledge_type = cell_id_to_ledge_type[cell_id];
    const int cell_redge_type = cell_id_to_redge_type[cell_id];

    int real_ledge_type = will_be_rotated == true ? cell_redge_type : cell_ledge_type;
    int real_redge_type = will_be_rotated == true ? cell_ledge_type : cell_redge_type;

    int num_total_edge_spacing_vio = 0;
    if(real_ledge_type > 0)
    {
      int spacing_between_redge_type1 = edge_spacing_rules[real_ledge_type + 1];
      int spacing_between_redge_type2 = edge_spacing_rules[real_ledge_type + 2];
  
      int num_vio_between_redge_type1 = 0;
      int num_vio_between_redge_type2 = 0;

      // Check spacing between right edge type 1
      for(int x = cur_x - spacing_between_redge_type1; x < cur_x; x++)
        for(int y = cur_y; y < cur_y + cell_h; y++)
          num_vio_between_redge_type1 += x >= 0 ? pixel_num_redge_type1[x * y_size + y] : 0;
  
      // Check spacing between right edge type 2
      for(int x = cur_x - spacing_between_redge_type2; x < cur_x; x++)
        for(int y = cur_y; y < cur_y + cell_h; y++)
          num_vio_between_redge_type2 += x >= 0 ? pixel_num_redge_type2[x * y_size + y] : 0;

      num_total_edge_spacing_vio += num_vio_between_redge_type1;
      num_total_edge_spacing_vio += num_vio_between_redge_type2;
    }
    
    if(real_redge_type > 0)
    {
      int spacing_between_ledge_type1 = edge_spacing_rules[real_redge_type + 1];
      int spacing_between_ledge_type2 = edge_spacing_rules[real_redge_type + 2];
  
      int num_vio_between_ledge_type1 = 0;
      int num_vio_between_ledge_type2 = 0;
  
      // Check spacing between left edge type 1
      for(int x = cur_x + cell_w; x < cur_x + cell_w + spacing_between_ledge_type1; x++)
        for(int y = cur_y; y < cur_y + cell_h; y++)
          num_vio_between_ledge_type1 += x < x_size ? pixel_num_ledge_type1[x * y_size + y] : 0;
  
      // Check spacing between left edge type 2
      for(int x = cur_x + cell_w; x < cur_x + cell_w + spacing_between_ledge_type2; x++)
        for(int y = cur_y; y < cur_y + cell_h; y++)
          num_vio_between_ledge_type2 += x < x_size ? pixel_num_ledge_type2[x * y_size + y] : 0;

      num_total_edge_spacing_vio += num_vio_between_ledge_type1;
      num_total_edge_spacing_vio += num_vio_between_ledge_type2;
    }

    cell_id_to_num_edge_spacing_vio[cell_id] = num_total_edge_spacing_vio;
  }
}

__global__ void countRowOrientVioKernel(
  const int  num_cells,
  const int* cell_ly_in_grid,
  const int* cell_h_in_grid,
  const int* cell_id_to_has_ground_at_bottom,
  const int* cell_need_rotation,
  const int* grid_y_to_is_vdd_up,
        int* cell_id_to_if_row_orient_violation)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int cur_y = cell_ly_in_grid[cell_id];
    const int cell_h = cell_h_in_grid[cell_id];
    const int cell_has_ground_at_bottom = cell_id_to_has_ground_at_bottom[cell_id];
    const int cell_rotated = cell_need_rotation[cell_id];
    const int is_row_vdd_up = grid_y_to_is_vdd_up[cur_y];

    int if_violation = 0;
    if((cell_h & 1) == 1)
    {
      if(cell_rotated == 0 and (cell_has_ground_at_bottom != is_row_vdd_up))
        if_violation = 1;
      else if(cell_rotated == 1 and (cell_has_ground_at_bottom == is_row_vdd_up))
        if_violation = 1;
    }
    else
    {
      if(cell_has_ground_at_bottom != is_row_vdd_up)
        if_violation = 1;
    }
    cell_id_to_if_row_orient_violation[cell_id] = if_violation;

    if(if_violation == 1)
      printf("RowViolation -> CellID: %d Y: %d\n", cell_id, cur_y);
  }
}

__device__ inline int runAlongAxis(
  const int  direction, 
  const int  cur_x, 
  const int  cur_y, 
  const int  x_size,
  const int  y_size,
  const int  group_id,
  const int  max_search,
  const int* pixel_id_to_valid,
  const int* pixel_id_to_group_id)
{
  int counter = 0;

  auto invalidity_test = [&] (int _x, int _y)
  {
    if(_x >= x_size or _x < 0 or _y >= y_size or _y < 0)
      return true;
    if(pixel_id_to_valid[_x * y_size + _y] != 1)
      return true;
    if(pixel_id_to_group_id[_x * y_size + _y] != group_id)
      return true;
    return false;
  };

  // 0: LEFT / 1: RIGHT / 2: UP / 3: DOWN
  if(direction == 0)
  {
    for(int x = cur_x; x >= cur_x - max_search; x--)
    {
      if(invalidity_test(x, cur_y) == true)
        break;
      counter++;
    }
  }
  else if(direction == 1)
  {
    for(int x = cur_x; x <= cur_x + max_search; x++)
    {
      if(invalidity_test(x, cur_y) == true)
        break;
      counter++;
    }
  }
  else if(direction == 2)
  {
    for(int y = cur_y; y <= cur_y + max_search; y++)
    {
      if(invalidity_test(cur_x, y) == true)
        break;
      counter++;
    }
  }
  else if(direction == 3)
  {
    for(int y = cur_y; y >= cur_y - max_search; y--)
    {
      if(invalidity_test(cur_x, y) == true)
        break;
      counter++;
    }
  }

  return counter;
}

__global__ void diagnosePixelKernel(
  const int   num_pixels,
  const int   x_size,
  const int   y_size,
  const int   max_search,
  const int   k_ratio,
  const int*  pixel_id_to_valid,
  const int*  pixel_id_to_group_id,
        int*  pixel_id_to_shape)
{
  const int pixel_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(pixel_id >= num_pixels)
    return;

  const int pixel_x = pixel_id / y_size;
  const int pixel_y = pixel_id % y_size;

  const int pixel_group_id = pixel_id_to_group_id[pixel_id];

  int left_counter 
    = runAlongAxis(0, pixel_x, pixel_y, x_size, y_size, 
        pixel_group_id, max_search, pixel_id_to_valid, pixel_id_to_group_id);

  int right_counter 
    = runAlongAxis(1, pixel_x, pixel_y, x_size, y_size, 
        pixel_group_id, max_search, pixel_id_to_valid, pixel_id_to_group_id);

  int up_counter 
    = runAlongAxis(2, pixel_x, pixel_y, x_size, y_size, 
        pixel_group_id, max_search, pixel_id_to_valid, pixel_id_to_group_id);

  int down_counter 
    = runAlongAxis(3, pixel_x, pixel_y, x_size, y_size, 
        pixel_group_id, max_search, pixel_id_to_valid, pixel_id_to_group_id);

  int local_width = left_counter + right_counter + 1;
  int local_height = up_counter + down_counter + 1;

  bool is_thin_h = local_width > k_ratio * local_height;
  bool is_thin_v = local_height > k_ratio * local_width;

  pixel_id_to_shape[pixel_id] = 
    pixel_id_to_valid[pixel_id] == 0 ? k_blockage :
    is_thin_h == true                ? k_shape_thin_h : 
    is_thin_v == true                ? k_shape_thin_v : k_shape_regular;
}

}
