#include <cstdio>
#include <memory>
#include <chrono>
#include <cassert>

#include "NesterovSolver.h"
#include "problem_instance/TargetFunction.h"

#include "Util.h"
#include "object/GPObject.h"
#include "cuda_linalg/CudaVectorAlgebra.h"

namespace skyplace 
{

__global__ void moveBackwardKernel(
  const int    num_var,
  const float  step_length, 
  const float* cur_pre_var,
  const float* cur_pre_grad,
        float* prv_pre_var)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < num_var)
    prv_pre_var[i] = cur_pre_var[i] - step_length * cur_pre_grad[i]; // move Backward
}

__global__ void moveForwardKernel(
  const int    num_var,
  const float  step_length, 
  const float  step_length_for_prediction, 
  const float* cur_var,
  const float* cur_pre_var,
  const float* cur_pre_grad,
        float* new_var,
        float* new_pre_var)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < num_var)
  {
    float next_x = cur_pre_var[i] + step_length * cur_pre_grad[i];
    float next_pre_x = next_x + step_length_for_prediction * (next_x - cur_var[i]);

    new_var[i] = next_x;
    new_pre_var[i] = next_pre_x;
  }
}

NesterovSolver::NesterovSolver(
  std::shared_ptr<HyperParam> param, 
  std::shared_ptr<TargetFunction> problem)
  : SolverBase(param, problem)
{
  curA_ = 0.0;
  step_length_ = 0.0;
  type_ = SolverType::NESTEROV;
}

void
NesterovSolver::initSolver()
{
  param_->printHyperParameters();

  initForCUDAKernel();
  
  // Step #1. Compute Initial Gradient
  target_function_->getInitialGrad(d_cur_pre_var_, d_cur_pre_grad_);

  // Step #2. Compute Initial Previous Coordinates
  // Since we don't have the previous predicted coordinates
  // at the first iteration, we have to make it up.
  moveBackward(param_->initOptCoef, d_cur_pre_var_, d_cur_pre_grad_, d_prv_pre_var_);
  // Compute the first Previous Predicted X / Y by going backward

  // Step #3. Update Pin Coordinates based on
  // Initial Previous Predicted Coordinates.
  // Compute the Previous Predicted Gradient 
  // based on the previous predicted Coordinate
  target_function_->updatePointAndGetGrad(d_prv_pre_var_, d_prv_pre_grad_);

  // Step #4. Compute the initial Step Length
  step_length_ = computeLipschitz(d_prv_pre_var_, 
                                  d_prv_pre_grad_,
                                  d_cur_pre_var_, 
                                  d_cur_pre_grad_);

  target_function_->updateParameters();
  //printf("Initial Nesterov StepLength: %.1E\n", step_length_);
}

void
NesterovSolver::solve()
{
  initSolver();

  printf("Nesterov Solve Start\n");
  target_function_->solveBgnCbk();

  int iter = 0;
  int backTrackIter = 0;
  curA_ = 1.0;

  auto solve_start_chrono = getChronoNow();
  for(; iter < param_->maxOptIter; iter++)
  {
    float prevA = curA_;
    curA_ = 0.5 * (1.0 + std::sqrt(4.0 * prevA * prevA + 1.0));
    float coeff = (prevA - 1.0) / curA_;

    step_length_ = backTracking(iter, coeff, backTrackIter);

    updateOneIteration(iter);

    if(target_function_->checkConvergence() == true)
      break;

    target_function_->iterEndCbk(iter, evalTime(solve_start_chrono), d_cur_var_);
  }

  target_function_->solveEndCbk(iter, evalTime(solve_start_chrono), d_cur_var_);
}

float
NesterovSolver::backTracking(int iter, float coeff, int& backTrackIter)
{
  backTrackIter = 0;
  float new_step_length = param_->minStepLength;

  for(; backTrackIter < param_->maxBackTrackIter; backTrackIter++)
  {
    moveForward(step_length_, 
                coeff,
                d_cur_var_,
                d_cur_pre_var_,
                d_cur_pre_grad_,
                d_new_var_,
                d_new_pre_var_);
  
    target_function_->updatePointAndGetGrad(d_new_pre_var_, d_new_pre_grad_);

    new_step_length = computeLipschitz(d_cur_pre_var_,
                                     d_cur_pre_grad_,
                                     d_new_pre_var_,
                                     d_new_pre_grad_);
  
    if(new_step_length >= step_length_ * 0.95)
      break;
    else if(new_step_length < param_->minStepLength) // minStepLength = 0.1
    {
      new_step_length = param_->minStepLength;
      break;
    }
  }
  return new_step_length;
}

void
NesterovSolver::moveForward(
  const float step_length, 
  const float step_length_for_prediction,
  const CudaVector<float>& d_cur_var,
  const CudaVector<float>& d_cur_pre_var,
  const CudaVector<float>& d_cur_direction,
        CudaVector<float>& d_new_var,
        CudaVector<float>& d_new_pre_var)
{
  int numThread = 64;
  int numBlockCell = (num_var_ - 1 + numThread) / numThread;

  moveForwardKernel<<<numBlockCell, numThread>>>(
    num_var_, 
    step_length, 
    step_length_for_prediction,
    d_cur_var.data(), 
    d_cur_pre_var.data(), 
    d_cur_direction.data(),
    d_new_var.data(), 
    d_new_pre_var.data());

  target_function_->clipToChipBoundary(d_new_var);
  target_function_->clipToChipBoundary(d_new_pre_var);
}

void
NesterovSolver::moveBackward(
  const float step_length, 
  const CudaVector<float>& d_cur_var,
  const CudaVector<float>& d_cur_direction,
        CudaVector<float>& d_new_var)
{
  int numThread = 64;
  int numBlockCell = (num_var_ - 1 + numThread) / numThread;

  moveBackwardKernel<<<numBlockCell, numThread>>>(
    num_var_, 
    step_length, 
    d_cur_var.data(), 
    d_cur_direction.data(),
    d_new_var.data());

  target_function_->clipToChipBoundary(d_new_var);
}

float
NesterovSolver::computeLipschitz(
  const CudaVector<float>& d_cur_pre_var, 
  const CudaVector<float>& d_cur_pre_grad,
  const CudaVector<float>& d_new_pre_var, 
  const CudaVector<float>& d_new_pre_grad)
{
  float delta_norm2          = computeDelta2Norm(d_cur_pre_var, d_new_pre_var, d_workSpaceForStepLength_);
  float coordi_distance      = delta_norm2 * delta_norm2;
  float coordi_distance_sqrt = std::sqrt(coordi_distance / static_cast<float>(num_var_));

  float grad_delta_norm2     = computeDelta2Norm(d_cur_pre_grad, d_new_pre_grad, d_workSpaceForStepLength_);
  float grad_distance        = grad_delta_norm2 * grad_delta_norm2;
  float grad_distance_sqrt   = std::sqrt(grad_distance / static_cast<float>(num_var_));

  return coordi_distance_sqrt / grad_distance_sqrt;
}

void
NesterovSolver::initForCUDAKernel()
{
  // Step #1. Vectors for Nesterov 
  d_workSpaceForStepLength_.resize(num_var_);

  d_prv_pre_var_.resize(num_var_);
  d_prv_pre_grad_.resize(num_var_);

  d_cur_pre_var_.resize(num_var_);
  d_cur_pre_grad_.resize(num_var_);

  d_new_pre_var_.resize(num_var_);
  d_new_pre_grad_.resize(num_var_);

  // Step #2. Set Initial Solution
  setInitialSolution();
}

void
NesterovSolver::setInitialSolution()
{
  // Host -> Device
  target_function_->exportToSolver(d_cur_var_);
  target_function_->exportToSolver(d_cur_pre_var_);
  // We don't have to initialize Previous X / Y
}

void
NesterovSolver::updateOneIteration(int iter)
{
  // Previous <= Current
  d_prv_pre_var_.swap(d_cur_pre_var_);
  d_prv_pre_grad_.swap(d_cur_pre_grad_);

  // Current <= Next
  d_cur_pre_var_.swap(d_new_pre_var_);
  d_cur_pre_grad_.swap(d_new_pre_grad_);

  d_cur_var_.swap(d_new_var_);
}

} // namespace skyplace
