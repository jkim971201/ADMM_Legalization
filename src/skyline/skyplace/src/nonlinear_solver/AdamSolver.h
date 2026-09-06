#ifndef ADAM_SOLVER_H
#define ADAM_SOLVER_H

#include "SolverBase.h"

namespace skyplace
{

class AdamSolver : public SolverBase
{
  public:

    AdamSolver();
    AdamSolver(
      std::shared_ptr<HyperParam> param, 
      std::shared_ptr<TargetFunction> gp_problem);

    // Main Loop
    void solve() override;

  private:

    void initSolver() override;
    void initForCUDAKernel() override;
    void setInitialSolution() override;
    void updateOneIteration(int iter) override;

    void moveForward(
      const float step_length,
      const CudaVector<float>& d_cur_var,
      const CudaVector<float>& d_cur_direction,
            CudaVector<float>& d_new_var);

    CudaVector<float> d_cur_direction_;
    CudaVector<float> d_cur_grad_;

    // Adam Optimizer
    float alpha_;
    float beta1_;
    float beta2_;

    float epsilon_;

    float beta1k_; // beta1_ ^ (k+1)
    float beta2k_; // beta2_ ^ (k+1)

    // First Moment
    CudaVector<float> d_cur_1st_momentum_;
    CudaVector<float> d_new_1st_momentum_;

    // Second Moment 
    CudaVector<float> d_cur_2nd_momentum_;
    CudaVector<float> d_new_2nd_momentum_;

    // Bias Corrected
    CudaVector<float> d_bias_corrected_1st_momentum_;
    CudaVector<float> d_bias_corrected_2nd_momentum_;

    void updateMoment(
      const CudaVector<float>& d_cur_1st_momentum,
      const CudaVector<float>& d_cur_2nd_momentum,
      const CudaVector<float>& d_cur_grad,
            CudaVector<float>& d_new_1st_momentum,
            CudaVector<float>& d_new_2nd_momentum);

    void updateDirection(
      const CudaVector<float>& d_new_1st_momentum,
      const CudaVector<float>& d_new_2nd_momentum,
            CudaVector<float>& d_bias_corrected_1st_momentum,
            CudaVector<float>& d_bias_corrected_2nd_momentum,
            CudaVector<float>& d_cur_direction);
};

}; // namespace skyline

#endif
