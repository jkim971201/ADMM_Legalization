#ifndef HYPERPARAM_H
#define HYPERPARAM_H

#include <cstdio>
#include <string>

namespace skyplace
{

enum SolverType
{
  NESTEROV,
  ADAM
};

enum InitialPlaceMethod
{
  RANDOM_START,
  SDP_RELAXATION
};

class HyperParam
{
  public:
    HyperParam()
    {
      maxOptIter       = 2500;
      maxBackTrackIter = 10;
      initLambda       = 1e-4;
      initGammaInvCoef = 0.01; 
      targetOverflow   = 0.07;
      minPhiCoef       = 0.95;
      maxPhiCoef       = 1.05;
      minPrecond       = 1.0;
      initOptCoef      = 100;
      referenceHpwl    = 200000;
      minStepLength    = 1.0;
      adam_alpha       = 100.0;
      adam_beta1       = 0.90;
      adam_beta2       = 0.999;
      log_freq         = 50;
      plot_mode        = false;
      plot_path        = "";
      plot_freq        = 5;
      solver_type      = NESTEROV;
      init_method      = RANDOM_START;
    }

    int maxOptIter;
    int maxBackTrackIter;

    float initLambda;
    float initGammaInvCoef;
    float targetOverflow;
    float minPhiCoef;
    float maxPhiCoef;
    float minPrecond;
    float initOptCoef;
    float referenceHpwl;
    float minStepLength;

    // Adam-related
    float adam_alpha;
    float adam_beta1;
    float adam_beta2;

    // Util Parameter
    int log_freq;
    bool plot_mode;
    int plot_freq;
    std::string plot_path;
    InitialPlaceMethod init_method;
    SolverType solver_type;

    void printHyperParameters() const
    {
      printf("\n");
      printf("Hyper Parameters\n");
      printf("  maxOptIter    : %12d   \n", maxOptIter);
      printf("  maxBackTrack  : %12d   \n", maxBackTrackIter);
      printf("  initLambda    : %12.6f \n", initLambda);
      printf("  initGammaInv  : %12.6f \n", initGammaInvCoef);
      printf("  TargetOvf     : %12.6f \n", targetOverflow);
      printf("  minPhiCoef    : %12.6f \n", minPhiCoef);
      printf("  maxPhiCoef    : %12.6f \n", maxPhiCoef);
      printf("  minPrecond    : %12.1f \n", minPrecond);
      printf("  initOptCoef   : %12.1f \n", initOptCoef);
      printf("  refHpwl       : %12.1f \n", referenceHpwl);
      printf("  minStepLength : %12.1f \n", minStepLength);
      printf("\n");
    }
};

}; // namespace skyplace

#endif 
