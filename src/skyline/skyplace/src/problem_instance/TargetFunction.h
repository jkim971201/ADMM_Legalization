#ifndef TARGET_FUNCTION_H
#define TARGET_FUNCTION_H

#include <vector>
#include <memory>

#include "cuda_linalg/CudaVector.h"
#include "painter/Painter.h" // No pimpl for unique_ptr

namespace skyplace 
{

class HyperParam;
class SkyPlaceDB;
class WireLengthGradient;
class DensityGradient;

using namespace cuda_linalg;

class TargetFunction
{
  public:

    // Constructor
    TargetFunction() {}
    TargetFunction(
      std::shared_ptr<SkyPlaceDB> db,
      std::shared_ptr<HyperParam> param);

    // APIs
    // var : {x_vector, y_vector, ... }
    
    void importFromDb();

    void updatePointAndGetGrad(
      const CudaVector<float>& var,
            CudaVector<float>& grad);

    void getInitialGrad(
      const CudaVector<float>& initial_var,
            CudaVector<float>& initial_grad);

    void clipToChipBoundary(CudaVector<float>& cell_pos);

    void exportToSolver(CudaVector<float>& var_from_solver); 
    // This will be used only to get initial solution

    void updateParameters();

    void solveBgnCbk();
    void solveEndCbk(int iter, double runtime, const CudaVector<float>& var);

    void iterBgnCbk(int iter);
    void iterEndCbk(int iter, double runtime, const CudaVector<float>& var);

    bool checkConvergence() const;

    int getNumVariable() const;

  private:

    int num_var_;
    int num_movable_;

    float gammaInv_;
    float lambda_;
    float minPrecond_;

    float hpwl_;
    float prevHpwl_;
    float overflow_;

    float x_min_;
    float y_min_;
    float x_max_;
    float y_max_;

    std::vector<float> h_cell_pos_;

    CudaVector<float> d_wl_grad_;
    CudaVector<float> d_density_grad_;
    CudaVector<float> d_cell2num_pin_;

    CudaVector<float> d_cell_real_width_;
    CudaVector<float> d_cell_real_height_;

    std::shared_ptr<SkyPlaceDB>         db_;
    std::shared_ptr<HyperParam>         param_;
    std::shared_ptr<WireLengthGradient> wireLength_;
    std::shared_ptr<DensityGradient>    density_;

    std::unique_ptr<Painter>            painter_;

    void prepareData();

    void updateGammaInv(float overflow);
    void updateLambda(float prevHpwl, float curHPWL);

    void diverge();

    void exportToDb(const CudaVector<float>& cell_pos);

    // Logging
    void printMetricRow() const;
    void printProgress(int iter, double runtime) const;
};

}; // namespace skyplace

#endif
