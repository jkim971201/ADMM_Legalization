#ifndef DENSITY_GRAD_H
#define DENSITY_GRAD_H

#include "cuda_linalg/CudaVector.h"

#include "SkyPlaceDB.h"
#include "PoissonSolver.h"

namespace skyplace
{

// C- Style Struct Version For CUDA Kernel Functions
typedef struct OVBIN
{
  int lxID = 0;
  int lyID = 0;
  int uxID = 0;
  int uyID = 0;
} OVBIN;

using namespace cuda_linalg;

class DensityGradient
{
  public:
    
    DensityGradient();
    DensityGradient(std::shared_ptr<SkyPlaceDB> db);

    // API
    void computeGrad(float* densityGradX,
                     float* densityGradY,
                     const float* cellCx,
                     const float* cellCy);

    // Getters
    int    numBinX        () const { return numBinX_;       }
    int    numBinY        () const { return numBinY_;       }
    float  overflow       () const { return overflow_;      }
    float  macroOverflow  () const { return macroOverflow_; }
    double densityTime    () const { return densityTime_;   }
    double poissonTime    () const { return poissonTime_;   }
    double binDenUpTime   () const { return binDenUpTime_;  }

    // To Visualize Density Information
    // These will be passed to Painter
    const float* getDevicePotential     () const { return d_ptr_binPotential_;          }
    const float* getDeviceBinDensity    () const { return d_ptr_binDensity_;            }
    const float* getDevicePreconditioner() const { return d_ptr_densityPreconditioner_; }

  private:

    int numBinX_;
    int numBinY_;
    int numMovable_;
    float sumMovableArea_;

    float overflow_;
    float macroOverflow_;
    float sumPenalty_;

    float dieLx_;
    float dieLy_;
    float dieUx_;
    float dieUy_;

    float targetDensity_;

    float binWidth_;
    float binHeight_;

    bool localLambdaMode_;

    std::shared_ptr<SkyPlaceDB> db_;
    std::unique_ptr<PoissonSolver> poissonSolver_;

    // Device Data
    CudaVector<float> d_cellDensityWidth_; 
    CudaVector<float> d_cellDensityHeight_; 
    CudaVector<float> d_cellDensityScale_;

    float* d_ptr_cellDensityWidth_;
    float* d_ptr_cellDensityHeight_;
    float* d_ptr_cellDensityScale_;

    CudaVector<int> d_isFiller_;
    int* d_ptr_isFiller_;

    CudaVector<int> d_isMacro_;
    int* d_ptr_isMacro_;

    CudaVector<float> d_fixedArea_;
    float* d_ptr_fixedArea_;

    CudaVector<float> d_macroArea_;
    float* d_ptr_macroArea_;

    CudaVector<float> d_scaledBinArea_;
    float* d_ptr_scaledBinArea_;

    CudaVector<float> d_binLx_;
    float* d_ptr_binLx_;

    CudaVector<float> d_binLy_;
    float* d_ptr_binLy_;

    CudaVector<float> d_binUx_;
    float* d_ptr_binUx_;

    CudaVector<float> d_binUy_;
    float* d_ptr_binUy_;

    CudaVector<float> d_movableArea_;
    float* d_ptr_movableArea_;

    CudaVector<float> d_fillerArea_;
    float* d_ptr_fillerArea_;

    CudaVector<float> d_overflowArea_;
    float* d_ptr_overflowArea_;

    CudaVector<float> d_macroOverflowArea_;
    float* d_ptr_macroOverflowArea_;

    CudaVector<float> d_binDensity_;
    float* d_ptr_binDensity_;

    CudaVector<float> d_binPenalty_;
    float* d_ptr_binPenalty_;

    CudaVector<float> d_binLambda_;
    float* d_ptr_binLambda_;

    CudaVector<float> d_densityPreconditioner_;
    float* d_ptr_densityPreconditioner_;

    CudaVector<float> d_binPotential_;
    float* d_ptr_binPotential_;

    CudaVector<float> d_electroForceX_;
    float* d_ptr_electroForceX_;

    CudaVector<float> d_electroForceY_;
    float* d_ptr_electroForceY_;

    // Runtime
    double densityTime_;
    double poissonTime_;
    double binDenUpTime_;

    void initForCUDAKernel();
};

} // namespace skyplace

#endif
