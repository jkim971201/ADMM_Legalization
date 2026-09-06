#ifndef QP_SOLVER_KERNELS_HPP
#define QP_SOLVER_KERNELS_HPP

#include "DeviceUtil.hpp"

namespace qp_solver
{

struct grid_to_dbu_float_functor
{
  int offset_;
  int unit_;

  grid_to_dbu_float_functor(int offset, int unit)
    : offset_(offset), unit_(unit) {}

  __host__ __device__
  float operator()(int grid) const
  {
    return float(grid * unit_ + offset_);
  }
};

struct slack_update_functor 
{
  float rho_;

  slack_update_functor(float rho): rho_(rho) {}

  __host__ __device__ 
  void operator()(thrust::tuple<float, float, float, float&> t) const
  { 
    float primal_product = thrust::get<0>(t);
    float bound = thrust::get<1>(t);
    float dual = thrust::get<2>(t);
    thrust::get<3>(t) = min(0.0f, - primal_product + bound - dual / rho_);
  }
};

struct dual_update_functor 
{
  float rho_;

  dual_update_functor(float rho): rho_(rho) {}

  __host__ __device__ 
  void operator()(thrust::tuple<float, float, float, float&> t) const
  { 
    float primal_product = thrust::get<0>(t);
    float slack = thrust::get<1>(t);
    float bound = thrust::get<2>(t);
    float dual = thrust::get<3>(t);
    thrust::get<3>(t) = dual + rho_ * (primal_product + slack - bound);
  }
};

struct absolute_diff
{
  __host__ __device__ 
  float operator()(thrust::tuple<float, float> t) const
  { 
    float x1 = thrust::get<0>(t);
    float x2 = thrust::get<1>(t);
    return fabsf(x1 - x2); 
  }
};

template<typename T>
__device__ inline void thomasAlgorithm(
  const int n, 
  const T   a, 
  const T   b, 
  const T   c_in, 
  const T*  d_in, 
  const T*  c_prime, 
        T*  d_prime, 
        T*  x)
{
  for(int i = 0; i < n; i++)
  {
    if(i == 0)
      d_prime[i] = d_in[i] / b;
    else
      d_prime[i] = (d_in[i] - a * d_prime[i - 1]) / (b - a * c_prime[i - 1]);
  }

  for(int i = 0; i < n; i++)
  {
    if(i == 0)
      x[n - i - 1] = d_prime[n - i - 1];
    else
      x[n - i - 1] = d_prime[n - i - 1] - c_prime[n - i - 1] * x[n - i];
  }
}

__global__ void assignLocalIndexKernel(
  const int    num_movable,
  const int    core_lx,
  const int    site_width,
  const int*   global_id_to_lx_in_grid,
  const int*   global_id_grouped_by_cut_row,
  const int*   global_id_to_original_lx,
        int*   local_id_to_global_id,
        float* local_id_to_primal,
        float* local_id_to_original_lx)
{
  const int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(thread_id >= num_movable)
    return;

  int local_id = thread_id;
  int global_id = global_id_grouped_by_cut_row[local_id];
  int original_lx_dbu = global_id_to_original_lx[global_id];

  local_id_to_global_id[local_id] = global_id;
  // convert to dbu (float)
  local_id_to_primal[local_id] = float(core_lx + site_width * global_id_to_lx_in_grid[global_id]);
  local_id_to_original_lx[local_id] = float(original_lx_dbu);
}

__global__ void generateConstraintsKernel(
  const int    num_cut_row,
  const int    core_lx,
  const int    site_width,
  const int*   cut_row_id_to_x_bgn,
  const int*   cut_row_id_to_x_end,
  const int*   cut_row_id_to_num_cells,
  const int*   cut_row_id_to_cell_offset,
  const int*   cell_grouped_by_cut_row,
  const int*   global_id_to_cell_width_in_grid,
  const int*   local_id_to_global_id,
        int*   local_id_to_pos_constraint_id,
        int*   local_id_to_neg_constraint_id,
        int*   constraint_id_to_pos_local_cell_id,
        int*   constraint_id_to_neg_local_cell_id,
        float* constraint_id_to_bound)
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id >= num_cut_row)
    return;

  const int num_cells = cut_row_id_to_num_cells[cut_row_id];
  const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];

  const int grid_x_bgn = cut_row_id_to_x_bgn[cut_row_id];
  const int grid_x_end = cut_row_id_to_x_end[cut_row_id];

  const int row_lx_dbu = site_width * grid_x_bgn + core_lx;
  const int row_ux_dbu = site_width * (grid_x_end + 1) + core_lx;

  for(int i = 0; i < num_cells + 1; i++)
  {
    int constraint_id = cell_offset + cut_row_id + i;
    int local_id  = (i == num_cells) ? -1 : i + cell_offset;
    int prev_cell_local_id = (i == 0) ? -1 : i + cell_offset - 1;
    int prev_cell_global_id = (i == 0) ? -1 : local_id_to_global_id[prev_cell_local_id];
    int prev_cell_width_dbu = (i == 0) ? 0 : site_width * global_id_to_cell_width_in_grid[prev_cell_global_id];

    if(local_id != -1)
      local_id_to_pos_constraint_id[local_id] = constraint_id;

    if(prev_cell_local_id != -1)
      local_id_to_neg_constraint_id[prev_cell_local_id] = constraint_id;

    constraint_id_to_pos_local_cell_id[constraint_id] = local_id;
    constraint_id_to_neg_local_cell_id[constraint_id] = prev_cell_local_id;

    if(num_cells == 0) // When cut_row has no cell
    {
      constraint_id_to_bound[constraint_id] = 0.0f;
    }
    else
    {
      constraint_id_to_bound[constraint_id] =
        (i == 0)         ? row_lx_dbu :
        (i == num_cells) ? prev_cell_width_dbu - row_ux_dbu : prev_cell_width_dbu;
    }
  }
}

__global__ void precomputeCPrimeKernel(
  const int    num_cut_row,
  const float  rho,
  const int*   cut_row_id_to_cell_offset,
  const int*   cut_row_id_to_num_cell,
        float* precomputed_c_prime)
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id >= num_cut_row)
    return;

  const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];
  const int num_cells_this_row = cut_row_id_to_num_cell[cut_row_id];
  for(int i = 0; i < num_cells_this_row; i++)
  {
    precomputed_c_prime[i + cell_offset] 
      = (i == 0) ? -rho / (1 + 2 * rho) : -rho / (1 + 2 * rho + rho * precomputed_c_prime[i + cell_offset - 1]);
  }
}

__global__ void updatePrimalXKernel(
  const int     num_cell,
  const float   rho, 
  const int*    cell_id_to_pos_const_id,
  const int*    cell_id_to_neg_const_id,
  const int*    cut_row_id_to_cell_offset,
  const int*    cut_row_id_to_num_cell,
  const float*  local_id_to_original_lx,
  const float*  constraint_id_to_bound,
  const float*  constraint_id_to_slack,
  const float*  constraint_id_to_dual,
        float*  thomas_workspace_d_in) /* return vector */
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cell)
  {
    float new_x = local_id_to_original_lx[cell_id];
    int const_id1 = cell_id_to_pos_const_id[cell_id];
    int const_id2 = cell_id_to_neg_const_id[cell_id];

    // compute [ x_init - rho * B^T * (y - b + lambda / rho) ]
    new_x -= rho * (constraint_id_to_slack[const_id1] - constraint_id_to_bound[const_id1] + constraint_id_to_dual[const_id1] / rho);
    new_x += rho * (constraint_id_to_slack[const_id2] - constraint_id_to_bound[const_id2] + constraint_id_to_dual[const_id2] / rho);
    thomas_workspace_d_in[cell_id] = new_x;
  }
}

__global__ void runThomasKernel(
  const int    num_row,
  const float  rho, 
  const int*   cut_row_id_to_cell_offset,
  const int*   cut_row_id_to_num_cell,
  const float* thomas_workspace_c_prime,
        float* thomas_workspace_d_in,
        float* thomas_workspace_d_prime,
        float* vector_x) /* return vector */
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id < num_row)
  {
    int cell_offset = cut_row_id_to_cell_offset[cut_row_id];
    int num_cell_this_row = cut_row_id_to_num_cell[cut_row_id];
  
    const float* c_prime = thomas_workspace_c_prime + cell_offset; // this is precomputed
    float* d_array = thomas_workspace_d_in    + cell_offset;
    float* d_prime = thomas_workspace_d_prime + cell_offset;
    float* x_array = vector_x + cell_offset;

    thomasAlgorithm(num_cell_this_row,   /* n       */
                    -rho,                /* a       */
                    1 + 2 * rho,         /* b       */
                    -rho,                /* c_in    */
                    d_array,             /* d_in    */
                    c_prime,             /* c_prime (workspace) */
                    d_prime,             /* d_prime (workspace) */
                    x_array);            /* solution x */
  }
}

__global__ void computeConstraintValueKernel(
  const int    num_constraint,
  const int*   constraint_id_to_pos_local_cell_id,
  const int*   constraint_id_to_neg_local_cell_id,
  const float* local_id_to_primal,
        float* constraint_id_to_primal_product) /* return vector */
{
  const int constraint_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(constraint_id >= num_constraint)
    return;

  const int pos_local_cell_id = constraint_id_to_pos_local_cell_id[constraint_id];
  const int neg_local_cell_id = constraint_id_to_neg_local_cell_id[constraint_id];

  if(pos_local_cell_id == -1 and neg_local_cell_id == -1)
    constraint_id_to_primal_product[constraint_id] = 0.0f;
  else if(pos_local_cell_id == -1 and neg_local_cell_id != -1)
    constraint_id_to_primal_product[constraint_id] = -local_id_to_primal[neg_local_cell_id];
  else if(pos_local_cell_id != -1 and neg_local_cell_id == -1)
    constraint_id_to_primal_product[constraint_id] = +local_id_to_primal[pos_local_cell_id];
  else
  {
    constraint_id_to_primal_product[constraint_id] 
      = +local_id_to_primal[pos_local_cell_id] - local_id_to_primal[neg_local_cell_id];
  }
}

__global__ void siteAlignNaiveKernel(
  const int    num_cut_row,
  const int    core_lx,
  const int    site_width,
  const int*   cut_row_id_to_x_bgn,
  const int*   cut_row_id_to_x_end,
  const int*   cut_row_id_to_num_cells,
  const int*   cut_row_id_to_cell_offset,
  const int*   global_id_to_cell_width_in_grid,
  const int*   local_id_to_global_id,
  const float* local_id_to_original_lx,
        float* local_id_to_primal,
        int*   global_id_to_lx_in_grid)
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id >= num_cut_row)
    return;

  const int num_cells = cut_row_id_to_num_cells[cut_row_id];
  const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];

  const int grid_x_bgn = cut_row_id_to_x_bgn[cut_row_id];
  const int grid_x_end = cut_row_id_to_x_end[cut_row_id];

  const int row_lx_dbu = site_width * grid_x_bgn + core_lx;
  const int row_ux_dbu = site_width * (grid_x_end + 1) + core_lx;

  int prev_cell_ux = row_lx_dbu;
  for(int i = 0; i < num_cells; i++)
  {
    const int local_id = i + cell_offset;
    const int global_id = local_id_to_global_id[local_id];
    const int width_dbu = global_id_to_cell_width_in_grid[global_id] * site_width;

    const int original_lx  = static_cast<int>(local_id_to_original_lx[local_id]);
    const int lx_clamp     = max(row_lx_dbu, static_cast<int>(local_id_to_primal[local_id]));
    const int dist         = lx_clamp - core_lx;
    const int new_lx_left  = (dist / site_width)     * site_width + core_lx;
    const int new_lx_right = (dist / site_width + 1) * site_width + core_lx;

    const float disp_left  = static_cast<float>(abs(new_lx_left - original_lx));
    const float disp_right = static_cast<float>(abs(new_lx_right - original_lx));

    int new_lx = new_lx_left;
    bool left_better = disp_left <= disp_right ? true : false;
    if(left_better == false)
    {
      if(i < num_cells - 1)
      {
        const int next_cell_lx = static_cast<int>(local_id_to_primal[local_id + 1]);
        const int next_cell_dist = next_cell_lx - core_lx;
        const int next_lx_left  = (next_cell_dist / site_width)     * site_width + core_lx;
        if(next_lx_left >= new_lx_right + width_dbu)
          new_lx = new_lx_right;
      }
      if(i == num_cells - 1)
      {
        if(new_lx_right + width_dbu <= row_ux_dbu)
          new_lx = new_lx_right;
      }
    }

    if(prev_cell_ux > new_lx)
      new_lx = new_lx_right;

    prev_cell_ux = new_lx + width_dbu;
    global_id_to_lx_in_grid[global_id] = (new_lx - core_lx) / site_width;
    local_id_to_primal[local_id] = static_cast<float>(new_lx);
  }
}

__global__ void siteAlignDPKernel(
    const int    num_cut_row,
    const int    core_lx,
    const int    site_width,
    const int*   cut_row_id_to_x_bgn,
    const int*   cut_row_id_to_x_end,
    const int*   cut_row_id_to_num_cells,
    const int*   cut_row_id_to_cell_offset,
    const int*   global_id_to_cell_width_in_grid,
    const int*   local_id_to_global_id,
    const float* local_id_to_original_lx,
          float* local_id_to_primal,
          int*   best_prev_state_when_move_left,  
          int*   best_prev_state_when_move_right, 
          int*   global_id_to_lx_in_grid)
{
  const int cut_row_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cut_row_id >= num_cut_row)
    return;

  const int num_cells = cut_row_id_to_num_cells[cut_row_id];
  if(num_cells == 0)
    return;

  const int row_lx_dbu = site_width * cut_row_id_to_x_bgn[cut_row_id] + core_lx;
  const int row_ux_dbu = site_width * (cut_row_id_to_x_end[cut_row_id] + 1) + core_lx;

  const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];

  constexpr int K_STATE_LEFT = 0;
  constexpr int K_STATE_RIGHT = 1;
  constexpr float K_OVERLAP_PENALTY = legalizer::k_infinity;

  float cost_when_prev_left = 0;
  float cost_when_prev_right = 0; 
  int ux_when_prev_left = row_lx_dbu;
  int ux_when_prev_right = row_lx_dbu; 
  for(int i = 0; i < num_cells; i++)
  {
    const int local_id     = i + cell_offset;
    const int global_id    = local_id_to_global_id[local_id];
    const int width_dbu    = global_id_to_cell_width_in_grid[global_id] * site_width;
    const int original_lx  = static_cast<int>(local_id_to_original_lx[local_id]);
    const int lx_clamp     = max(row_lx_dbu, static_cast<int>(local_id_to_primal[local_id]));
    const int dist         = lx_clamp - core_lx;
    const int new_lx_left  = (dist / site_width)     * site_width + core_lx;
    const int new_lx_right = (dist / site_width + 1) * site_width + core_lx;
    const float disp_left  = static_cast<float>(abs(new_lx_left - original_lx));
    const float disp_right = static_cast<float>(abs(new_lx_right - original_lx));

    float cost_when_move_left; 
    float cost_when_move_right; 
    if(i == 0)
    {
      float ov_left = (new_lx_left < row_lx_dbu) ? K_OVERLAP_PENALTY : 0;
      float ov_right = (new_lx_right < row_lx_dbu or new_lx_right + width_dbu > row_ux_dbu) ? K_OVERLAP_PENALTY : 0;
      cost_when_move_left = disp_left + ov_left;
      cost_when_move_right = disp_right + ov_right;

      // no meaning since no backtracking when i = 0
      best_prev_state_when_move_left[local_id] = -1; 
      best_prev_state_when_move_right[local_id] = -1;
    }
    else
    {
      // STATE_LEFT
      float ov_when_move_left_prev_left = (new_lx_left < ux_when_prev_left) ? K_OVERLAP_PENALTY : 0;
      float ov_when_move_left_prev_right = (new_lx_left < ux_when_prev_right) ? K_OVERLAP_PENALTY : 0;
      float cost_when_move_left_prev_left = cost_when_prev_left + disp_left + ov_when_move_left_prev_left;
      float cost_when_move_left_prev_right = cost_when_prev_right + disp_left + ov_when_move_left_prev_right;
      if(cost_when_move_left_prev_left <= cost_when_move_left_prev_right)
      {
        cost_when_move_left = cost_when_move_left_prev_left; 
        best_prev_state_when_move_left[local_id] = K_STATE_LEFT; 
      }
      else                  
      {
        cost_when_move_left = cost_when_move_left_prev_right;
        best_prev_state_when_move_left[local_id] = K_STATE_RIGHT; 
      }

      // STATE_RIGHT
      float bound_penalty = (new_lx_right + width_dbu > row_ux_dbu) ? K_OVERLAP_PENALTY : 0;
      float ov_when_move_right_prev_left = ((new_lx_right < ux_when_prev_left) ? K_OVERLAP_PENALTY : 0) + bound_penalty;
      float ov_when_move_right_prev_right = ((new_lx_right < ux_when_prev_right) ? K_OVERLAP_PENALTY : 0) + bound_penalty;
      float cost_when_move_right_prev_left = cost_when_prev_left + disp_right + ov_when_move_right_prev_left;
      float cost_when_move_right_prev_right = cost_when_prev_right + disp_right + ov_when_move_right_prev_right;
      if(cost_when_move_right_prev_left <= cost_when_move_right_prev_right)
      {
        cost_when_move_right = cost_when_move_right_prev_left; 
        best_prev_state_when_move_right[local_id] = K_STATE_LEFT; 
      }
      else                      
      {
        cost_when_move_right = cost_when_move_right_prev_right; 
        best_prev_state_when_move_right[local_id] = K_STATE_RIGHT; 
      }
    }
    
    cost_when_prev_left = cost_when_move_left;
    cost_when_prev_right = cost_when_move_right;
    ux_when_prev_left = new_lx_left + width_dbu;
    ux_when_prev_right = new_lx_right + width_dbu;
  }

  // Backtracking
  int best_state = (cost_when_prev_left <= cost_when_prev_right) ? K_STATE_LEFT : K_STATE_RIGHT;
  for(int i = num_cells - 1; i >= 0; i--)
  {
    const int local_id  = i + cell_offset;
    const int global_id = local_id_to_global_id[local_id];
    const int lx_clamp  = max(row_lx_dbu, static_cast<int>(local_id_to_primal[local_id]));
    const int dist      = lx_clamp - core_lx;
    global_id_to_lx_in_grid[global_id] = (best_state == K_STATE_LEFT) ? dist / site_width : dist / site_width + 1;
    int prev_state = (best_state == K_STATE_LEFT) ? best_prev_state_when_move_left[local_id] : best_prev_state_when_move_right[local_id];
    best_state = prev_state;
  }
}

__global__ void debugKernel(
  const int    num_cut_row,
  const int    core_lx,
  const int    site_width,
  const int*   cut_row_id_to_num_cells,
  const int*   cut_row_id_to_cell_offset,
  const int*   local_id_to_pos_constraint_id,
  const int*   local_id_to_neg_constraint_id,
  const float* local_id_to_primal,
  const int*   constraint_id_to_pos_local_cell_id,
  const int*   constraint_id_to_neg_local_cell_id,
  const float* constraint_id_to_bound,
  const float* constraint_id_to_primal_product)
{
  for(int i = 0; i < num_cut_row; i++)
  {
    const int cut_row_id = i;
    const int num_cells = cut_row_id_to_num_cells[cut_row_id];
    const int cell_offset = cut_row_id_to_cell_offset[cut_row_id];

    printf("cut_row_id: %d\n", cut_row_id);
    for(int j = 0; j < num_cells + 1; j++)
    {
      const int constraint_id = cell_offset + cut_row_id + j;
      const int local_id = (j == num_cells) ? -1 : j + cell_offset;
      const int prev_local_id = (j == 0) ? -1 : j + cell_offset - 1;

      const float cell_lx_dbu = (local_id == num_cells) ? 0 : local_id_to_primal[local_id];
      const float prev_cell_lx_dbu = (local_id == 0) ? 0 : local_id_to_primal[prev_local_id];

      const float bound = constraint_id_to_bound[constraint_id];
      const float primal_product = constraint_id_to_primal_product[constraint_id];
      const float residual = primal_product - bound;

      //if(cut_row_id == 3)
      {
        printf("local_id: %d prev_local_id: %d cell_lx: %f prev_cell_lx: %f primal_product: %f bound: %f residual: %f\n", 
          local_id, prev_local_id, cell_lx_dbu, prev_cell_lx_dbu, primal_product, bound, residual);

        if(residual < 0.0f)
          printf("NegResidual!! (%f)\n", residual);
      }
    }
  }
}

}

#endif
