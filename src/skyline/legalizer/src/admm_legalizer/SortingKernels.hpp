#ifndef ADMM_SORTING_KERNEL_HPP
#define ADMM_SORTING_KERNEL_HPP

#include "GlobalUtil.h"

namespace admm_legalizer
{

__global__ void computeKeyArrayForSortingKernel(
  const int  x_size,
  const int  y_size,
  const int  num_cells,
  const int  coeff_random,
  const int  coeff_ovf,
  const int  coeff_area,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* pixel_usage,
  const int* partition_id_to_num_cell,
  const int* partition_id_to_offset,
  const int* cell_random_priority,
  const int* cell_id_grouped_by_partition,
        int* cell_key_grouped_by_partition)
{
  const int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_cells)
    return;

  const int id_in_partition = thread_id;
  const int cell_id = cell_id_grouped_by_partition[thread_id];

  const int w_in_grid = cell_w_in_grid[cell_id];
  const int h_in_grid = cell_h_in_grid[cell_id];
  const int lx_in_grid = cell_lx_in_grid[cell_id];
  const int ly_in_grid = cell_ly_in_grid[cell_id];

  const int x1 = lx_in_grid;
  const int y1 = ly_in_grid;
  const int x2 = lx_in_grid + w_in_grid;
  const int y2 = ly_in_grid + h_in_grid;

  int sum_ovf = 0;
  for(int x = x1; x < x2; x++)
    for(int y = y1; y < y2; y++)
      sum_ovf += max(pixel_usage[x * y_size + y] - 1, 0);

  const int random_key = cell_random_priority[cell_id];
  const int ovf_key = sum_ovf;
  const int area_key = w_in_grid * h_in_grid;

  cell_key_grouped_by_partition[id_in_partition] 
    = coeff_random * random_key + coeff_ovf * ovf_key + coeff_area * area_key;
}

__global__ void sortEachPartitionKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int  buffer_size, // each partition must have cells less than this.
  const int  coeff_random,
  const int  coeff_ovf,
  const int  coeff_area,
  const int* partition_id_to_offset,
  const int* pixel_usage,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_random_priority,
        int* cell_id_grouped_by_partition)
{
  // shared_mem size is 2 * buffer_size
  extern __shared__ int merged_buffer[]; 

  int* buffer_key = merged_buffer;
  int* buffer_val = merged_buffer + buffer_size;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int partition_offset = partition_id_to_offset[partition_id];
  const int num_cells_this_partition = partition_id_to_offset[partition_id + 1] - partition_offset;

  int* cells_this_partition = cell_id_grouped_by_partition + partition_offset;

  assert(num_cells_this_partition <= buffer_size);

  // Initialize the buffer as INT_MIN
  for(int i = num_cells_this_partition + thread_id; i < buffer_size; i += stride)
  {
    buffer_key[i] = k_int_min;
    buffer_val[i] = -1;
  }

  __syncthreads();

  // First compute overflow sum
  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    const int cell_id = cells_this_partition[i];
    const int w_in_grid = cell_w_in_grid[cell_id];
    const int h_in_grid = cell_h_in_grid[cell_id];
    const int lx_in_grid = cell_lx_in_grid[cell_id];
    const int ly_in_grid = cell_ly_in_grid[cell_id];

    const int x1 = lx_in_grid;
    const int y1 = ly_in_grid;
    const int x2 = lx_in_grid + w_in_grid;
    const int y2 = ly_in_grid + h_in_grid;

    int sum_ovf = 0;
    for(int x = x1; x < x2; x++)
      for(int y = y1; y < y2; y++)
        sum_ovf += max(pixel_usage[x * y_size + y] - 1, 0);

    const int random_key = cell_random_priority[cell_id];
    const int ovf_key = sum_ovf;
    const int area_key = w_in_grid * h_in_grid;

    buffer_key[i] = coeff_random * random_key + coeff_ovf * ovf_key + coeff_area * area_key;
    buffer_val[i] = cell_id;

    assert(cell_id < num_cells);
  }

  __syncthreads();

  // Sort by key (overflow)
  bitonicSortDescending(stride, thread_id, buffer_size, buffer_key, buffer_val);

  __syncthreads();

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    assert(buffer_val[i] < num_cells);
    cells_this_partition[i] = buffer_val[i];
  }
}

__device__ inline void mergeSubroutine(
  const int  left, 
  const int  mid, 
  const int  right,
  const int* key,
  const int* val,
        int* sorted_key,
        int* sorted_val)
{
  int i = left;
  int j = mid;
  int k = left;

  while(i < mid && j < right) 
  {
    if(key[i] > key[j])
    {
      sorted_key[k] = key[i];
      sorted_val[k] = val[i];
      i++;
    }
    else
    {
      sorted_key[k] = key[j];
      sorted_val[k] = val[j];
      j++;
    }
    k++;
  }

  while(i < mid) 
  {
    sorted_key[k] = key[i];
    sorted_val[k] = val[i];
    k++;
    i++;
  }

  while(j < right)  
  {
    sorted_key[k] = key[j];
    sorted_val[k] = val[j];
    k++;
    j++;
  }
}

__device__ inline void mergeSortDescending(
  int stride, int thread_id, int n, int* key_sorted, int* val_sorted, int* buffer_key, int* buffer_val)
{
  for(int w = 1; w < n; w *= 2) 
  {
    for(int j = thread_id; j < n / 2; j += stride)
    {
      int left  = j * 2 * w;
      int mid   = min(left + w, n);
      int right = min(left + 2 * w, n);

      mergeSubroutine(
        left, mid, right,
        key_sorted, 
        val_sorted, 
        buffer_key, 
        buffer_val);
    }

    __syncthreads(); 

    for(int i = thread_id; i < n; i += stride)
    {
      key_sorted[i] = buffer_key[i];
      val_sorted[i] = buffer_val[i];
    }

    __syncthreads();
  }
}

__global__ void bitonicSortKernel(
  const int  buffer_size,
  const int* partition_id_to_offset,
        int* cell_id_grouped_by_partition,
        int* cell_key_grouped_by_partition)
{
  // shared_mem size is 2 * buffer_size
  extern __shared__ int merged_buffer[]; 

  int* buffer_key = merged_buffer;
  int* buffer_val = merged_buffer + buffer_size;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int partition_offset = partition_id_to_offset[partition_id];
  const int num_cells_this_partition = partition_id_to_offset[partition_id + 1] - partition_offset;

  int* key_this_partition = cell_key_grouped_by_partition + partition_offset;
  int* cells_this_partition = cell_id_grouped_by_partition + partition_offset;

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    buffer_key[i] = key_this_partition[i];
    buffer_val[i] = cells_this_partition[i];
  }

  for(int i = num_cells_this_partition + thread_id; i < buffer_size; i += stride)
  {
    buffer_key[i] = k_int_min;
    buffer_val[i] = -1;
  }

  __syncthreads();

  // Sort by key (overflow)
  bitonicSortDescending(stride, thread_id, buffer_size, buffer_key, buffer_val);

  __syncthreads();

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    key_this_partition[i] = buffer_key[i];
    cells_this_partition[i] = buffer_val[i];
  }
}

__global__ void mergeSortKernel(
  const int  buffer_size,
  const int* partition_id_to_offset,
        int* cell_id_grouped_by_partition,
        int* cell_key_grouped_by_partition)
{
  // shared_mem size is 2 * buffer_size
  extern __shared__ int merged_buffer[]; 

  int* buffer_key = merged_buffer;
  int* buffer_val = merged_buffer + buffer_size;

  int* buffer_key_sorted = merged_buffer + 2 * buffer_size;
  int* buffer_val_sorted = merged_buffer + 3 * buffer_size;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int partition_offset = partition_id_to_offset[partition_id];
  const int num_cells_this_partition = partition_id_to_offset[partition_id + 1] - partition_offset;

  int* key_this_partition = cell_key_grouped_by_partition + partition_offset;
  int* cells_this_partition = cell_id_grouped_by_partition + partition_offset;

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    buffer_key[i] = key_this_partition[i];
    buffer_val[i] = cells_this_partition[i];

    buffer_key_sorted[i] = buffer_key[i];
    buffer_val_sorted[i] = buffer_val[i];
  }

  for(int i = num_cells_this_partition + thread_id; i < buffer_size; i += stride)
  {
    buffer_key[i] = k_int_min;
    buffer_val[i] = -1;

    buffer_key_sorted[i] = buffer_key[i];
    buffer_val_sorted[i] = buffer_val[i];
  }

  __syncthreads();

  // Sort by key (overflow)
  mergeSortDescending(
    stride, 
    thread_id, 
    buffer_size, 
    buffer_key_sorted, 
    buffer_val_sorted, 
    buffer_key, 
    buffer_val);

  __syncthreads();

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    key_this_partition[i] = buffer_key[i];
    cells_this_partition[i] = buffer_val[i];
  }
}

__global__ void sortEachPartitionWithMergeSortKernel(
  const int  num_cells,
  const int  x_size,
  const int  y_size,
  const int  buffer_size, // each partition must have cells less than this.
  const int  coeff_random,
  const int  coeff_ovf,
  const int  coeff_area,
  const int* partition_id_to_offset,
  const int* pixel_usage,
  const int* cell_w_in_grid,
  const int* cell_h_in_grid,
  const int* cell_lx_in_grid,
  const int* cell_ly_in_grid,
  const int* cell_random_priority,
        int* cell_id_grouped_by_partition)
{
  // shared_mem size is 2 * buffer_size
  extern __shared__ int merged_buffer[]; 

  int* buffer_key = merged_buffer;
  int* buffer_val = merged_buffer + buffer_size;

  int* buffer_key_sorted = merged_buffer + 2 * buffer_size;
  int* buffer_val_sorted = merged_buffer + 3 * buffer_size;

  const int partition_id = blockIdx.x;
  const int stride = blockDim.x;
  const int thread_id = threadIdx.x;

  const int partition_offset = partition_id_to_offset[partition_id];
  const int num_cells_this_partition = partition_id_to_offset[partition_id + 1] - partition_offset;

  int* cells_this_partition = cell_id_grouped_by_partition + partition_offset;

  assert(num_cells_this_partition <= buffer_size);

  // Initialize the buffer as INT_MIN
  for(int i = num_cells_this_partition + thread_id; i < buffer_size; i += stride)
  {
    buffer_key[i] = k_int_min;
    buffer_val[i] = -1;

    buffer_key_sorted[i] = k_int_min;
    buffer_val_sorted[i] = -1;
  }

  __syncthreads();

  // First compute overflow sum
  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    const int cell_id = cells_this_partition[i];
    const int w_in_grid = cell_w_in_grid[cell_id];
    const int h_in_grid = cell_h_in_grid[cell_id];
    const int lx_in_grid = cell_lx_in_grid[cell_id];
    const int ly_in_grid = cell_ly_in_grid[cell_id];

    const int x1 = lx_in_grid;
    const int y1 = ly_in_grid;
    const int x2 = lx_in_grid + w_in_grid;
    const int y2 = ly_in_grid + h_in_grid;

    int sum_ovf = 0;
    for(int x = x1; x < x2; x++)
      for(int y = y1; y < y2; y++)
        sum_ovf += max(pixel_usage[x * y_size + y] - 1, 0);

    const int random_key = cell_random_priority[cell_id];
    const int ovf_key = sum_ovf;
    const int area_key = w_in_grid * h_in_grid;

    buffer_key[i] = coeff_random * random_key + coeff_ovf * ovf_key + coeff_area * area_key;
    buffer_val[i] = cell_id;

    buffer_key_sorted[i] = buffer_key[i];
    buffer_val_sorted[i] = buffer_val[i];

    assert(cell_id < num_cells);
  }

  __syncthreads();

  // Sort by key (overflow)
  mergeSortDescending(
    stride, 
    thread_id, 
    buffer_size, 
    buffer_key_sorted, 
    buffer_val_sorted, 
    buffer_key, 
    buffer_val);

  __syncthreads();

  for(int i = thread_id; i < num_cells_this_partition; i += stride)
  {
    assert(buffer_val[i] < num_cells);
    cells_this_partition[i] = buffer_val_sorted[i];
  }
}

__global__ void debugCubSortKernel(
  const int  num_partition,
  const int* partition_id_to_offset,
  const int* cell_key_grouped_by_partition,
  const int* cell_id_grouped_by_partition)
{
  for(int partition_id = 0; partition_id < num_partition; partition_id++)
  {
    const int partition_offset = partition_id_to_offset[partition_id];
    const int num_cells_this_partition = partition_id_to_offset[partition_id + 1] - partition_offset;

    const int* key_this_partition = cell_key_grouped_by_partition + partition_offset;
    const int* cells_this_partition = cell_id_grouped_by_partition + partition_offset;

    printf("Key\n");
    for(int i = 0; i < num_cells_this_partition; i++)
      printf("%d ", key_this_partition[i]);
    printf("\n");

    printf("Value\n");
    for(int i = 0; i < num_cells_this_partition; i++)
      printf("%d ", cells_this_partition[i]);
    printf("\n");
  }
}

}

#endif
