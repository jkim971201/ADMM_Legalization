#ifndef NESTEROV_SOLVER_H
#define NESTEROV_SOLVER_H

#include "SolverBase.h"

namespace skyplace
{

class NesterovSolver : public SolverBase
{
  public:

    NesterovSolver();
    NesterovSolver(
      std::shared_ptr<HyperParam> param, 
      std::shared_ptr<TargetFunction> problem);

    // Main Loop
    void solve() override;

  private:

    void initSolver() override;
    void initForCUDAKernel() override;
    void setInitialSolution() override;
    void updateOneIteration(int iter) override;

    float backTracking(int iter, float coeff, int& backTrackIter);

    // Self-Adaptive Parameters for Nesterov Acceleration
    float curA_;
    float step_length_;

    // Previous Predicted Coordinates
    CudaVector<float> d_prv_pre_var_;
    CudaVector<float> d_prv_pre_grad_;

    // Current Predicted Coordinates
    CudaVector<float> d_cur_pre_var_;
    CudaVector<float> d_cur_pre_grad_;

    // Next Predicted Coordinates
    CudaVector<float> d_new_pre_var_;
    CudaVector<float> d_new_pre_grad_;

    CudaVector<float> d_workSpaceForStepLength_;

    // GPU-related Functions
    void moveForward(
      const float step_length, 
      const float step_length_for_prediction,
      const CudaVector<float>& d_cur_var,
      const CudaVector<float>& d_cur_pre_var,
      const CudaVector<float>& d_cur_direction,
            CudaVector<float>& d_new_var,
            CudaVector<float>& d_new_pre_var);

    void moveBackward(
      const float step_length, 
      const CudaVector<float>& d_cur_var,
      const CudaVector<float>& d_cur_direction,
            CudaVector<float>& d_new_var);

    float computeLipschitz(
      const CudaVector<float>& d_cur_var,
      const CudaVector<float>& d_cur_pre_grad,
      const CudaVector<float>& d_new_pre_var, 
      const CudaVector<float>& d_new_pre_grad);
};

}; // namespace skyplace 

#endif
