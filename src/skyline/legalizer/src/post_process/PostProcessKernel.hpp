#include "DeviceUtil.hpp"

namespace legalizer
{

__global__ void rotateKernel(
  const int  num_cells,
  const int* grid_y_to_is_vdd_up,
  const int* cell_ly_in_grid,
  const int* cell_h_in_grid,
  const int* cell_has_ground_at_bottom,
        int* cell_need_rotation)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int cell_h = cell_h_in_grid[cell_id];
    const int cur_y = cell_ly_in_grid[cell_id];
    const int is_row_vdd_up = grid_y_to_is_vdd_up[cur_y];
    const int cell_vss_bottom = cell_has_ground_at_bottom[cell_id];

    if((cell_h & 1) == 1)
      cell_need_rotation[cell_id] = (cell_vss_bottom == is_row_vdd_up) ? 0 : 1;
    // do nothing for even-heighted cells
  }
}

__global__ void stampCellsKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
        int* pixel_to_cell_id)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    const int x_min = cell_lx_in_grid[cell_id];
    const int y_min = cell_ly_in_grid[cell_id];
    const int x_max = x_min + cell_w_in_grid[cell_id];
    const int y_max = y_min + cell_h_in_grid[cell_id];
    for(int x = x_min; x < x_max; x++)
    {
      for(int y = y_min; y < y_max; y++)
      {
        assert(x >= 0 and x < x_size);
        assert(y >= 0 and y < y_size);
        pixel_to_cell_id[x * y_size + y] = cell_id;
      }
    }
  }
}

__global__ void countCutRowsKernel(
  const int  x_size,
  const int  y_size,
  const int* cell_h_in_grid,
  const int* pixel_id_to_group_id,
  const int* pixel_id_to_valid,
  const int* pixel_id_to_cell_id,
        int* grid_y_to_num_cut_rows)
{
  const int y = blockIdx.x * blockDim.x + threadIdx.x;
  if(y < y_size)
  {
    int prev_group_id = -999;
    int prev_is_valid = -999;
    int prev_is_single = -999;
    for(int x = 0; x < x_size; x++)
    {
      int pixel_id = x * y_size + y;
      int group_id = pixel_id_to_group_id[pixel_id];
      int is_valid = pixel_id_to_valid[pixel_id];
      int cell_id  = pixel_id_to_cell_id[pixel_id];
      int cell_h   = (cell_id >= 0) ? cell_h_in_grid[cell_id] : 0;
      int is_single = (cell_h <= 1) ? 1 : 0;

      // Check if the current segment ends
      bool segment_ends = (prev_group_id != group_id) || (prev_is_valid != is_valid) || (cell_h > 1);
      if(segment_ends == true)
        if(prev_is_valid == 1 and prev_is_single == 1)
          grid_y_to_num_cut_rows[y] += 1;
    
      // Finish the current segment in the rightmost pixel
      if(x == x_size - 1)
        if(segment_ends == false and is_valid == 1 and cell_h <= 1)
          grid_y_to_num_cut_rows[y] += 1;
    
      prev_group_id = group_id;
      prev_is_valid = is_valid;
      prev_is_single = is_single;
    }
  }
}

__global__ void generateCutRowsKernel(
  const int  x_size,
  const int  y_size,
  const int* cell_h_in_grid,
  const int* pixel_id_to_group_id,
  const int* pixel_id_to_valid,
  const int* pixel_id_to_cell_id,
  const int* grid_y_to_cut_row_offset,
        int* cut_row_id_to_grid_y,
        int* cut_row_id_to_grid_x_bgn,
        int* cut_row_id_to_grid_x_end)
{
  const int y = blockIdx.x * blockDim.x + threadIdx.x;
  if(y < y_size)
  {
    const int offset_this_y = grid_y_to_cut_row_offset[y];
    int num_already_made_in_this_y = 0;

    int prev_group_id = -999;
    int prev_is_valid = -999;
    int prev_is_single = -999;

    int grid_x_bgn = 0;
    int grid_x_end = 0;

    auto create_cut_row = [&] (int x_bgn, int x_end)
    {
      int cut_row_id = offset_this_y + num_already_made_in_this_y;
      cut_row_id_to_grid_y[cut_row_id] = y;
      cut_row_id_to_grid_x_bgn[cut_row_id] = x_bgn;
      cut_row_id_to_grid_x_end[cut_row_id] = x_end;
      num_already_made_in_this_y++;
    };

    for(int x = 0; x < x_size; x++)
    {
      int pixel_id = x * y_size + y;
      int group_id = pixel_id_to_group_id[pixel_id];
      int is_valid = pixel_id_to_valid[pixel_id];
      int cell_id  = pixel_id_to_cell_id[pixel_id];
      int cell_h   = (cell_id >= 0) ? cell_h_in_grid[cell_id] : 0;
      int is_single = (cell_h <= 1) ? 1 : 0;
      
      // Check if the current segment ends
      bool segment_ends = (prev_group_id != group_id) || (prev_is_valid != is_valid) || (prev_is_single != is_single);   
      if(segment_ends == true)
      {
        grid_x_end = x - 1;
        if(prev_is_valid == 1 and prev_is_single == 1)
          create_cut_row(grid_x_bgn, grid_x_end);
        grid_x_bgn = x;
      }

      // Finish the current segment in the rightmost pixel
      if(x == x_size - 1 and segment_ends == false and is_valid == 1 and cell_h <= 1)
        create_cut_row(grid_x_bgn, x_size - 1);

      prev_group_id = group_id;
      prev_is_valid = is_valid;
      prev_is_single = is_single;
    }
  }
}

__global__ void markCutRowToGridKernel(
  const int  x_size,
  const int  y_size,
  const int  num_cut_rows,
  const int* cut_row_id_to_grid_y,
  const int* cut_row_id_to_grid_x_bgn,
  const int* cut_row_id_to_grid_x_end,
        int* pixel_id_to_cut_row_id)
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id < num_cut_rows)
  {
    int cut_row_y = cut_row_id_to_grid_y[cut_row_id];
    int cut_row_x_bgn = cut_row_id_to_grid_x_bgn[cut_row_id];
    int cut_row_x_end = cut_row_id_to_grid_x_end[cut_row_id];
    for(int x = cut_row_x_bgn; x <= cut_row_x_end; x++)
      pixel_id_to_cut_row_id[x * y_size + cut_row_y] = cut_row_id;
  }
}

__global__ void markCutRowToCellKernel(
  const int  num_cells,
  const int  num_cut_rows,
  const int  x_size,
  const int  y_size,
  const int* pixel_id_to_cut_row_id,
  const int* cell_h_in_grid,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
        int* cell_id_to_cut_row_id)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cells)
  {
    int cell_h = cell_h_in_grid[cell_id];
    int cell_lx = cell_lx_in_grid[cell_id];
    int cell_ly = cell_ly_in_grid[cell_id];

    int pixel_id = cell_lx * y_size + cell_ly;
    int cut_row_id = pixel_id_to_cut_row_id[pixel_id];

    cell_id_to_cut_row_id[cell_id] = (cell_h == 1) ? cut_row_id : num_cut_rows;

    // Mark invalid cut_row_id to multi-height cells.
    if(cut_row_id == -1 and cell_h == 1)
    {
      printf("CUTROWID: %d CELLID: %d x: %d y: %d\n", cut_row_id, cell_id, cell_lx, cell_ly);
      assert(0);
    }
  }
}

__device__ inline void bitonicSortAscending(
  int stride, int thread_id, int n, int* buffer_key, int* buffer_value)
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
      swap(buffer_value, i, j);
    }
  };

  auto descent_swap = [&] (int i, int j)
  {
    if(buffer_key[i] < buffer_key[j])
    {
      swap(buffer_key, i, j);
      swap(buffer_value, i, j);
    }
  };

  for(int k = 2; k <= n; k *= 2)
  {
    for(int j = k / 2; j > 0; j /= 2) 
    {
      for(int i = thread_id; i < n; i += stride)
      {
        int ij = i ^ j;
        if(ij < i) 
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

__global__ void sortCellsEachCutRowKernel(
  const int  buffer_size,
  const int* cut_row_id_to_num_cells,
  const int* cut_row_id_to_cell_offset,
  const int* cell_id_to_grid_lx,
        int* cell_grouped_by_cut_row)
{
  const int cut_row_id = blockIdx.x;

  // shared_mem size is 2 * buffer_size (max_num_cell_in_cut_row)
  extern __shared__ int merged_buffer[]; 

  int* buffer_key = merged_buffer;
  int* buffer_value = merged_buffer + buffer_size;

  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];
  const int num_cells_this_cut_row = cut_row_id_to_num_cells[cut_row_id];

  int* cells_this_cut_row = cell_grouped_by_cut_row + cell_offset;
  assert(num_cells_this_cut_row <= buffer_size);

  for(int i = thread_id; i < buffer_size; i += stride)
  {
    // Initialize padded part as INT_MIN
    if(i >= num_cells_this_cut_row)
    {
      buffer_key[i] = k_int_max;
      buffer_value[i] = -1;
    }
    else
    {
      int cell_id = cells_this_cut_row[i];
      buffer_key[i] = cell_id_to_grid_lx[cell_id];
      buffer_value[i] = cell_id;
    }
  }

  __syncthreads();

  // Sort by key (overflow)
  bitonicSortAscending(stride, thread_id, buffer_size, buffer_key, buffer_value);

  __syncthreads();

  for(int i = thread_id; i < num_cells_this_cut_row; i += stride)
    cells_this_cut_row[i] = buffer_value[i];
}

__global__ void debugSortKernel(
  const int  num_cut_row,
  const int* cut_row_id_to_num_cells,
  const int* cut_row_id_to_cell_offset,
  const int* cut_row_id_to_grid_y,
  const int* cut_row_id_to_grid_x_bgn,
  const int* cut_row_id_to_grid_x_end,
  const int* cell_id_to_grid_lx,
  const int* cell_id_to_grid_ly,
  const int* cell_id_to_grid_w,
  const int* cell_grouped_by_cut_row)
{
  for(int cut_row_id = 0; cut_row_id < num_cut_row; cut_row_id++)
  {
    int num_cells = cut_row_id_to_num_cells[cut_row_id];
    int cell_offset = cut_row_id_to_cell_offset[cut_row_id];
    int cut_row_y = cut_row_id_to_grid_y[cut_row_id];
    int cut_row_x_bgn = cut_row_id_to_grid_x_bgn[cut_row_id];
    int cut_row_x_end = cut_row_id_to_grid_x_end[cut_row_id];
    printf("CutRowID: %5d x: %4d, %4d y: %5d\n", cut_row_id, cut_row_x_bgn, cut_row_x_end, cut_row_y);
    for(int i = 0; i < num_cells; i++)
    {
      int cell_id = cell_grouped_by_cut_row[i + cell_offset];
      int cell_lx = cell_id_to_grid_lx[cell_id];
      int cell_ly = cell_id_to_grid_ly[cell_id];
      int cell_w  = cell_id_to_grid_w[cell_id];
      printf("  CellID: %5d (%5d, %5d) Width: %2d\n", cell_id, cell_lx, cell_ly, cell_w);
    }
  }
}

}
