#include <cstdio>
#include <memory>
#include <chrono>
#include <cmath>

#include "AdamSolver.h"
#include "problem_instance/TargetFunction.h"

#include "Util.h"
#include "object/GPObject.h"
#include "cuda_linalg/CudaVectorAlgebra.h"

namespace skyplace
{

__global__ void moveForwardKernelAdam(
  const int    num_var,
  const float  step_length,
  const float* cur_var,
  const float* cur_direction,
        float* new_var)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < num_var)
    new_var[i] = cur_var[i] - step_length * cur_direction[i];
}

__global__ void updateDirectionKernelAdam(
  const int    num_var,
  const float  epsilon,
  const float* bias_corrected_1st_momentum,
  const float* bias_corrected_2nd_momentum,
        float* cur_direction)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < num_var)
  {
    float bc_1st = bias_corrected_1st_momentum[i];
    float bc_2nd = max(0.0, bias_corrected_2nd_momentum[i]);
    cur_direction[i] = bc_1st / ( sqrtf(bc_2nd) + epsilon );
  }
}

AdamSolver::AdamSolver() : alpha_(100.0), beta1_(0.9), beta2_(0.999) {}

AdamSolver::AdamSolver(
  std::shared_ptr<HyperParam>     param, 
  std::shared_ptr<TargetFunction> gp_problem)
  : SolverBase(param, gp_problem)
{
  alpha_  = param->adam_alpha;
  beta1_  = param->adam_beta1;
  beta2_  = param->adam_beta2;
  beta1k_ = beta1_;
  beta2k_ = beta2_;
  type_   = SolverType::ADAM;

  initForCUDAKernel();
}

void
AdamSolver::initSolver()
{
  // Step #1. Reset parameters
  alpha_  = param_->adam_alpha; // learning_rate
  beta1k_ = beta1_;
  beta2k_ = beta2_;

  // Step #2. Set Initial Solution
  setInitialSolution();

  target_function_->getInitialGrad(d_cur_var_, d_cur_grad_);

  target_function_->updateParameters();

  d_cur_1st_momentum_.fillZero();
  d_cur_2nd_momentum_.fillZero();
}

void
AdamSolver::moveForward(
  const float step_length,
  const CudaVector<float>& d_cur_var,
  const CudaVector<float>& d_cur_direction,
        CudaVector<float>& d_new_var)
{
  int num_thread = 64;
  int num_block_cell = (num_var_ - 1 + num_thread) / num_thread;

  moveForwardKernelAdam<<<num_block_cell, num_thread>>>(
    num_var_, 
    step_length,
    d_cur_var.data(), 
    d_cur_direction.data(),
    d_new_var.data() );

  target_function_->clipToChipBoundary(d_new_var);
}

void
AdamSolver::updateMoment(
  const CudaVector<float>& d_cur_1st_momentum,
  const CudaVector<float>& d_cur_2nd_momentum,
  const CudaVector<float>& d_cur_grad,
        CudaVector<float>& d_new_1st_momentum,
        CudaVector<float>& d_new_2nd_momentum)
{
  const float one_minus_beta1 = (1.0 - beta1_);
  const float one_minus_beta2 = (1.0 - beta2_);

  vectorAdd(beta1_, -one_minus_beta1, d_cur_1st_momentum, d_cur_grad, d_new_1st_momentum);
  vectorAxpySquare(beta2_, one_minus_beta2, d_cur_2nd_momentum, d_cur_grad, d_new_2nd_momentum);
}

void
AdamSolver::updateDirection(
  const CudaVector<float>& d_new_1st_momentum,
  const CudaVector<float>& d_new_2nd_momentum,
        CudaVector<float>& d_bias_corrected_1st_momentum,
        CudaVector<float>& d_bias_corrected_2nd_momentum,
        CudaVector<float>& d_cur_direction)
{
  float coeff1 = 1.0 / (1.0 - beta1k_);
  float coeff2 = 1.0 / (1.0 - beta2k_);

  vectorMulScalar(coeff1, d_new_1st_momentum, d_bias_corrected_1st_momentum);
  vectorMulScalar(coeff2, d_new_2nd_momentum, d_bias_corrected_2nd_momentum);

  int num_thread = 64;
  int num_block_cell = (num_var_ - 1 + num_thread) / num_thread;

  updateDirectionKernelAdam<<<num_block_cell, num_thread>>>(
    num_var_, 
    epsilon_,
    d_bias_corrected_1st_momentum.data(),
    d_bias_corrected_2nd_momentum.data(),
    d_cur_direction.data());
}

void
AdamSolver::solve()
{
  initSolver();

  printf("Adam Solve Start\n");
  target_function_->solveBgnCbk();

  auto solve_start_chrono = getChronoNow();

  int iter = 0;
  for(; iter < param_->maxOptIter; iter++)
  {
    target_function_->iterBgnCbk(iter);

    // Step #1: Compute gradient
    target_function_->updatePointAndGetGrad(d_cur_var_, d_cur_grad_);

    // Step #2: Update moment
    updateMoment(d_cur_1st_momentum_,
                 d_cur_2nd_momentum_,
                 d_cur_grad_,
                 d_new_1st_momentum_,
                 d_new_2nd_momentum_);

    // Step #3: Update direction
    updateDirection(d_new_1st_momentum_,
                    d_new_2nd_momentum_,
                    d_bias_corrected_1st_momentum_,
                    d_bias_corrected_2nd_momentum_,
                    d_cur_direction_);
    
    // Step #4: Move forward with the direction vector
    moveForward(alpha_, d_cur_var_, d_cur_direction_, d_new_var_);

    updateOneIteration(iter);

    if(target_function_->checkConvergence() == true)
      break;

    target_function_->iterEndCbk(iter, evalTime(solve_start_chrono), d_cur_var_);
  }

  target_function_->solveEndCbk(iter, evalTime(solve_start_chrono), d_cur_var_);
}

void
AdamSolver::initForCUDAKernel()
{
  // Step #1. Vectors for Adam
  d_cur_grad_.resize(num_var_);

  d_cur_direction_.resize(num_var_);

  d_cur_1st_momentum_.resize(num_var_);
  d_cur_2nd_momentum_.resize(num_var_);

  d_new_1st_momentum_.resize(num_var_);
  d_new_2nd_momentum_.resize(num_var_);

  d_bias_corrected_1st_momentum_.resize(num_var_);
  d_bias_corrected_2nd_momentum_.resize(num_var_);
}

void
AdamSolver::setInitialSolution()
{
  // Host -> Device
  target_function_->exportToSolver(d_cur_var_);
}

void
AdamSolver::updateOneIteration(int iter)
{
  // Current <= Next
  d_cur_var_.swap(d_new_var_);
  d_cur_1st_momentum_.swap(d_new_1st_momentum_);
  d_cur_2nd_momentum_.swap(d_new_2nd_momentum_);

  alpha_ *= 0.997;

  beta1k_ *= beta1_;
  beta2k_ *= beta2_;
}

}; // namespace skyplace
