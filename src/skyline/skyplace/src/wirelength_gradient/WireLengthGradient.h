#ifndef WIRELENGTH_GRAD_H
#define WIRELENGTH_GRAD_H

#include "cuda_linalg/CudaVector.h"
#include "SkyPlaceDB.h"

#include <limits>

namespace skyplace
{

using namespace cuda_linalg;

constexpr float k_float_max = std::numeric_limits<float>::max();

class WireLengthGradient
{
  public:
    
    WireLengthGradient();
    WireLengthGradient(std::shared_ptr<SkyPlaceDB> db);

    void computeGrad(float* wlGradX, float* wlGradY);

    void updatePinCoordinates(const float* d_cellCx, const float* d_cellCy);

    float computeHPWL();

		// Getters
    double wlGradTime() const { return wl_grad_time_; }

		// Setters
    void setGammaInv (const float gammaInv) { gammaInv_ = gammaInv; }

  private:

    std::shared_ptr<SkyPlaceDB> db_;

    int numNet_;
    int numPin_;
    int numCell_;

    float gammaInv_;

    CudaVector<float> d_pinX_;
    CudaVector<float> d_pinY_;

    float* d_ptr_pinX_;
    float* d_ptr_pinY_;
  
    CudaVector<int> d_pin2Net_;
    int* d_ptr_pin2Net_;

    CudaVector<int> d_pin2Cell_;
    int* d_ptr_pin2Cell_;
    
    CudaVector<int> d_netStart_;
    int* d_ptr_netStart_;

    CudaVector<float> d_pinOffsetX_;
    float* d_ptr_pinOffsetX_;

    CudaVector<float> d_pinOffsetY_;
    float* d_ptr_pinOffsetY_;

    CudaVector<float> d_netBBoxWidth_;
    float* d_ptr_netBBoxWidth_;

    CudaVector<float> d_netBBoxHeight_;
    float* d_ptr_netBBoxHeight_;

    CudaVector<float> d_waForEachNetX_;
    float* d_ptr_waForEachNetX_;

    CudaVector<float> d_waForEachNetY_;
    float* d_ptr_waForEachNetY_;

    CudaVector<float> d_maxPinX_;
    CudaVector<float> d_minPinX_;

    CudaVector<float> d_maxPinY_;
    CudaVector<float> d_minPinY_;

    float* d_ptr_maxPinX_;
    float* d_ptr_minPinX_;

    float* d_ptr_maxPinY_;
    float* d_ptr_minPinY_;

    CudaVector<float> d_apX_;
    CudaVector<float> d_bpX_;
    CudaVector<float> d_cpX_;

    CudaVector<float> d_apY_;
    CudaVector<float> d_bpY_;
    CudaVector<float> d_cpY_;

    CudaVector<float> d_amX_;
    CudaVector<float> d_bmX_;
    CudaVector<float> d_cmX_;

    CudaVector<float> d_amY_;
    CudaVector<float> d_bmY_;
    CudaVector<float> d_cmY_;

    float* d_ptr_apX_;
    float* d_ptr_bpX_;
    float* d_ptr_cpX_;

    float* d_ptr_apY_;
    float* d_ptr_bpY_;
    float* d_ptr_cpY_;

    float* d_ptr_amX_;
    float* d_ptr_bmX_;
    float* d_ptr_cmX_;

    float* d_ptr_amY_;
    float* d_ptr_bmY_;
    float* d_ptr_cmY_;

    CudaVector<float> d_pinGradX_;
    CudaVector<float> d_pinGradY_;

    float* d_ptr_pinGradX_;
    float* d_ptr_pinGradY_;

    CudaVector<float> d_netWeight_;
    float* d_ptr_netWeight_;

    double wl_grad_time_;

    // Methods
    void initForCUDAKernel();
};

}; // namespace skyplace 

#endif
