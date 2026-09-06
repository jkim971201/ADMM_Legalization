#include "LGHyperParameters.h"
#include "database/CudaDatabase.h"
#include "CudaQpSolver.h"
#include "QpSolverKernels.hpp"
#include "DeviceUtil.hpp"

#include "cuda_linalg/CudaVectorAlgebra.h"

#include <thrust/transform_reduce.h>

namespace qp_solver
{

CudaQpSolver::CudaQpSolver(
  const std::shared_ptr<LGHyperParameters> hyperparameter,
        std::shared_ptr<CudaDatabase> database)
  : database_(database), params_(hyperparameter)
{
  rho_ = params_->qp_admm_rho;
}

void
CudaQpSolver::initializeSolver(
  const cuda_linalg::CudaVector<int>& cut_row_id_to_grid_x_bgn,
  const cuda_linalg::CudaVector<int>& cut_row_id_to_grid_x_end,
  const cuda_linalg::CudaVector<int>& cut_row_id_to_cell_offset,
  const cuda_linalg::CudaVector<int>& cut_row_id_to_num_cells,
  const cuda_linalg::CudaVector<int>& cell_grouped_by_cut_row)
{
  num_cut_row_ = static_cast<int>(cut_row_id_to_num_cells.size());
  num_movable_ = computeVectorSum(cut_row_id_to_num_cells);
  num_constraint_ = num_movable_ + num_cut_row_;

  global_id_to_lx_in_grid_ = database_->getCellLxInGrid();

  local_id_to_global_id_.resize(num_movable_);
  local_id_to_pos_constraint_id_.resize(num_movable_);
  local_id_to_neg_constraint_id_.resize(num_movable_);
  local_id_to_primal_.resize(num_movable_);
  local_id_to_original_lx_.resize(num_movable_);

  cut_row_id_to_cell_offset_.resize(num_cut_row_);
  cut_row_id_to_num_cells_.resize(num_cut_row_);
  cut_row_id_to_grid_x_bgn_.resize(num_cut_row_);
  cut_row_id_to_grid_x_end_.resize(num_cut_row_);

  cut_row_id_to_cell_offset_ = cut_row_id_to_cell_offset;
  cut_row_id_to_num_cells_ = cut_row_id_to_num_cells;
  cut_row_id_to_grid_x_bgn_ = cut_row_id_to_grid_x_bgn;
  cut_row_id_to_grid_x_end_ = cut_row_id_to_grid_x_end;

  constraint_id_to_residual_.resize(num_constraint_);
  constraint_id_to_primal_product_.resize(num_constraint_);
  constraint_id_to_slack_.resize(num_constraint_);
  constraint_id_to_dual_.resize(num_constraint_);
  constraint_id_to_bound_.resize(num_constraint_);
  constraint_id_to_pos_local_cell_id_.resize(num_constraint_);
  constraint_id_to_neg_local_cell_id_.resize(num_constraint_);

  const int core_lx = database_->getCoreLx();
  const int site_width = database_->getSiteWidth();
  const auto& global_id_to_lx_in_grid = database_->getCellLxInGrid();
  const auto& global_id_to_original_lx = database_->getCellOriginalLxInDbu();

  const int num_thread = 128;
  const int num_block_cell = (num_movable_ + num_thread - 1) / num_thread;

  assignLocalIndexKernel<<<num_block_cell, num_thread>>>(
    num_movable_,
    core_lx,
    site_width,
    global_id_to_lx_in_grid.data(),
    cell_grouped_by_cut_row.data(),
    global_id_to_original_lx.data(),
    local_id_to_global_id_.data(),
    local_id_to_primal_.data(),
    local_id_to_original_lx_.data());

  const int num_block_row = (num_cut_row_ + num_thread - 1) / num_thread;

  const auto& global_id_to_cell_width_in_grid = database_->getCellWidthInGrid();

  generateConstraintsKernel<<<num_block_row, num_thread>>>(
    num_cut_row_,
    core_lx,
    site_width,
    cut_row_id_to_grid_x_bgn.data(),
    cut_row_id_to_grid_x_end.data(),
    cut_row_id_to_num_cells.data(),
    cut_row_id_to_cell_offset.data(),
    cell_grouped_by_cut_row.data(),
    global_id_to_cell_width_in_grid.data(),
    local_id_to_global_id_.data(),
    local_id_to_pos_constraint_id_.data(),
    local_id_to_neg_constraint_id_.data(),
    constraint_id_to_pos_local_cell_id_.data(),
    constraint_id_to_neg_local_cell_id_.data(),
    constraint_id_to_bound_.data());

  thomas_precomputed_c_prime_.resize(num_movable_);
  thomas_workspace_c_prime_.resize(num_movable_);
  thomas_workspace_d_prime_.resize(num_movable_);
  thomas_workspace_d_in_.resize(num_movable_);

  // Precompute c'
  precomputeCPrimeKernel<<<num_block_row, num_thread>>>(
    num_cut_row_,
    rho_,
    cut_row_id_to_cell_offset.data(),
    cut_row_id_to_num_cells.data(),
    thomas_precomputed_c_prime_.data());

  //printf("NumCutRow    : %d\n", num_cut_row_);
  //printf("NumMovable   : %d\n", num_movable_);
  //printf("NumConstraint: %d\n", num_constraint_);
}

bool
CudaQpSolver::checkTermination(int iter)
{
  vectorAdd(1.0f, -1.0f, 
    constraint_id_to_primal_product_, 
    constraint_id_to_bound_, 
    constraint_id_to_residual_);

  float min_residual = cuda_linalg::computeVectorMin(constraint_id_to_residual_);
  //printf("MinRes: %f\n", min_residual);
  bool is_convergence = (min_residual < -1 or iter < params_->qp_admm_min_iter) ? false : true;;
  return is_convergence;
}

void
CudaQpSolver::updatePrimalX()
{
  const int num_thread = 128;
  const int num_block1 = (num_movable_ + num_thread - 1) / num_thread;
  const int num_block2 = (num_cut_row_ + num_thread - 1) / num_thread;

  updatePrimalXKernel<<<num_block1, num_thread>>>(
    num_movable_,
    rho_,
    local_id_to_pos_constraint_id_.data(),
    local_id_to_neg_constraint_id_.data(),
    cut_row_id_to_cell_offset_.data(),
    cut_row_id_to_num_cells_.data(),
    local_id_to_original_lx_.data(),
    constraint_id_to_bound_.data(),
    constraint_id_to_slack_.data(),
    constraint_id_to_dual_.data(),
    thomas_workspace_d_in_.data());

  runThomasKernel<<<num_block2, num_thread>>>(
    num_cut_row_,
    rho_,
    cut_row_id_to_cell_offset_.data(),
    cut_row_id_to_num_cells_.data(),
    thomas_precomputed_c_prime_.data(),
    thomas_workspace_d_in_.data(),
    thomas_workspace_d_prime_.data(),
    local_id_to_primal_.data());
}

void
CudaQpSolver::updatePrimalY()
{
  auto zip_bgn 
    = thrust::make_zip_iterator(thrust::make_tuple(
        constraint_id_to_primal_product_.begin(),
        constraint_id_to_bound_.begin(),
        constraint_id_to_dual_.begin(),
        constraint_id_to_slack_.begin()));
        
  auto zip_end 
    = thrust::make_zip_iterator(thrust::make_tuple(
        constraint_id_to_primal_product_.end(),
        constraint_id_to_bound_.end(),
        constraint_id_to_dual_.end(),
        constraint_id_to_slack_.end()));

  thrust::for_each(zip_bgn, zip_end, slack_update_functor(rho_));
}

void
CudaQpSolver::updateDual()
{
  auto zip_bgn 
    = thrust::make_zip_iterator(thrust::make_tuple(
        constraint_id_to_primal_product_.begin(),
        constraint_id_to_slack_.begin(),
        constraint_id_to_bound_.begin(),
        constraint_id_to_dual_.begin()));

  auto zip_end
    = thrust::make_zip_iterator(thrust::make_tuple(
        constraint_id_to_primal_product_.end(),
        constraint_id_to_slack_.end(),
        constraint_id_to_bound_.end(),
        constraint_id_to_dual_.end()));

  thrust::for_each(zip_bgn, zip_end, dual_update_functor(rho_));
}

void
CudaQpSolver::computeConstraintValue()
{
  const int num_thread = 128;
  const int num_block  = (num_constraint_ + num_thread - 1) / num_thread;

  computeConstraintValueKernel<<<num_block, num_thread>>>(
    num_constraint_,
    constraint_id_to_pos_local_cell_id_.data(),
    constraint_id_to_neg_local_cell_id_.data(),
    local_id_to_primal_.data(),
    constraint_id_to_primal_product_.data());

  /*
  debugKernel<<<1, 1>>>(
    num_cut_row_,
    database_->getCoreLx(),
    database_->getSiteWidth(),
    cut_row_id_to_num_cells_.data(),
    cut_row_id_to_cell_offset_.data(),
    local_id_to_pos_constraint_id_.data(),
    local_id_to_neg_constraint_id_.data(),
    local_id_to_primal_.data(),
    constraint_id_to_pos_local_cell_id_.data(),
    constraint_id_to_neg_local_cell_id_.data(),
    constraint_id_to_bound_.data(),
    constraint_id_to_primal_product_.data());
    */
}

void 
CudaQpSolver::solve()
{
  //float site_width_float = static_cast<float>(database_->getSiteWidth());
  //float initial_hor_displace = computeDisplace(local_id_to_primal_) / site_width_float;
  //float initial_primal_objective = computePrimalObjective(local_id_to_primal_);
  //printf("Displace  before QP ADMM: %5.2f (sites)\n", initial_hor_displace);
  //printf("PrimalObj before QP ADMM: %f\n", initial_primal_objective);

  for(int admm_iter = 0; admm_iter < params_->qp_admm_max_iter; admm_iter++)
  {
    /* ===================== Main ADMM Iteration ====================== */
    /* 1. Solve subproblem w.r.t. x_k */
    updatePrimalX();

    computeConstraintValue();
    
    /* 2. Solve subproblem w.r.t. y_k */
    updatePrimalY();
   
    /* 3. Update Dual Variable lambda */
    updateDual();
  
    //float primal_obj = computePrimalObjective(local_id_to_primal_);
    //printf("Iter: %3d PrimalObjective: %f\n", admm_iter, primal_obj);

    bool is_termination = checkTermination(admm_iter);
    if(is_termination == true)
    {
      //float primal_obj = computePrimalObjective(local_id_to_primal_);
      printf("QP ADMM converged at iter %d\n", admm_iter);
      break;
    }
    /* ================================================================ */
  }

  if(params_->qp_admm_dp_site_align == true)
    siteAlignDP();
  else
    siteAlignNaive();

  //float final_hor_displace = computeDisplace(local_id_to_primal_) / site_width_float;
  //float final_primal_objective = computePrimalObjective(local_id_to_primal_);
  //printf("Displace   after QP ADMM: %5.2f (sites)\n", final_hor_displace);
  //printf("PrimalObj  after QP ADMM: %f\n", final_primal_objective);

  const auto& cell_ly_in_grid = database_->getCellLyInGrid();
  database_->commitPlaceInGrid(global_id_to_lx_in_grid_, cell_ly_in_grid);
}

float
CudaQpSolver::computePrimalObjective(const cuda_linalg::CudaVector<float>& primal) const
{
  float norm2 = cuda_linalg::computeDelta2NormSquare(primal, local_id_to_original_lx_);
  return norm2;
}

float
CudaQpSolver::computeDisplace(const cuda_linalg::CudaVector<float>& primal) const
{
  auto zip_bgn = thrust::make_zip_iterator(thrust::make_tuple(
    primal.begin(), local_id_to_original_lx_.begin()));

  auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(
    primal.end(), local_id_to_original_lx_.end()));

  float displace_hor_dbu 
    = thrust::transform_reduce(zip_bgn, zip_end, absolute_diff(), 0.0f, thrust::plus<float>());
  return displace_hor_dbu;
}

void
CudaQpSolver::siteAlignNaive()
{
  const int num_thread = 128;
  const int num_block  = (num_cut_row_ + num_thread - 1) / num_thread;

  const auto& global_id_to_cell_width_in_grid = database_->getCellWidthInGrid();

  siteAlignNaiveKernel<<<num_block, num_thread>>>(
    num_cut_row_,
    database_->getCoreLx(),
    database_->getSiteWidth(),
    cut_row_id_to_grid_x_bgn_.data(),
    cut_row_id_to_grid_x_end_.data(),
    cut_row_id_to_num_cells_.data(),
    cut_row_id_to_cell_offset_.data(),
    global_id_to_cell_width_in_grid.data(),
    local_id_to_global_id_.data(),
    local_id_to_original_lx_.data(),
    local_id_to_primal_.data(),
    global_id_to_lx_in_grid_.data());
}

void
CudaQpSolver::siteAlignDP()
{
  const int num_thread = 128;
  const int num_block  = (num_cut_row_ + num_thread - 1) / num_thread;

  const auto& global_id_to_cell_width_in_grid = database_->getCellWidthInGrid();

  cuda_linalg::CudaVector<int> best_prev_state_when_move_left(num_movable_);
  cuda_linalg::CudaVector<int> best_prev_state_when_move_right(num_movable_);

  siteAlignDPKernel<<<num_block, num_thread>>>(
    num_cut_row_,
    database_->getCoreLx(),
    database_->getSiteWidth(),
    cut_row_id_to_grid_x_bgn_.data(),
    cut_row_id_to_grid_x_end_.data(),
    cut_row_id_to_num_cells_.data(),
    cut_row_id_to_cell_offset_.data(),
    global_id_to_cell_width_in_grid.data(),
    local_id_to_global_id_.data(),
    local_id_to_original_lx_.data(),
    local_id_to_primal_.data(),
    best_prev_state_when_move_left.data(),
    best_prev_state_when_move_right.data(),
    global_id_to_lx_in_grid_.data());
}

}
