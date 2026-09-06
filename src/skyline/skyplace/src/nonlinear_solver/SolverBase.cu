#include "SolverBase.h"
#include "problem_instance/TargetFunction.h"

namespace skyplace
{

SolverBase::SolverBase() : num_var_(0) {}

SolverBase::SolverBase(
  std::shared_ptr<HyperParam>     param, 
  std::shared_ptr<TargetFunction> problem)
{
  param_ = param;
  target_function_ = problem;

  initializeSolverBase(target_function_);
}

void
SolverBase::initializeSolverBase(std::shared_ptr<TargetFunction> target_function)
{
  num_var_ = target_function->getNumVariable();
  d_cur_var_.resize(num_var_);
  d_new_var_.resize(num_var_);
}

SolverType
SolverBase::getSolverType() const
{
  return type_;
}

}
