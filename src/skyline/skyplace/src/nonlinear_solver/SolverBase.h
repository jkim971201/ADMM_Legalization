#ifndef SOLVER_BASE_H 
#define SOLVER_BASE_H

#include <memory>

#include "cuda_linalg/CudaVector.h"
#include "HyperParam.h"

namespace skyplace
{

using namespace cuda_linalg;

class TargetFunction;

class SolverBase 
{
  public:

    SolverBase();
    SolverBase(
      std::shared_ptr<HyperParam> param, 
      std::shared_ptr<TargetFunction> gp_problem);

    SolverType getSolverType() const;

    // Pure virtual function
    virtual void solve() = 0;

  protected:

    // Pure virtual function
    virtual void initSolver() = 0;
    virtual void initForCUDAKernel() = 0;
    virtual void setInitialSolution() = 0;
    virtual void updateOneIteration(int iter) = 0;

    void initializeSolverBase(std::shared_ptr<TargetFunction> target_function);

    int num_var_;

    SolverType type_;

		std::shared_ptr<HyperParam>     param_;
    std::shared_ptr<TargetFunction> target_function_;

    /* Device Array */
    CudaVector<float> d_cur_var_;
    CudaVector<float> d_new_var_;
};

} // namespace skyplace 

#endif
