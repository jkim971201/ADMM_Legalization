#ifndef SKYPLACE_H
#define SKYPLACE_H

#include <string>
#include <memory>

namespace db {
  class dbDatabase;
}

namespace skyplace 
{

class SkyPlaceDB;
class InitialPlacer;
class TargetFunction;
class DensityGradient;
class WireLengthGradient;
class SolverBase;
class HyperParam;

class SkyPlace 
{
  public:

    SkyPlace();
    SkyPlace(std::shared_ptr<db::dbDatabase> db);
 
    ~SkyPlace();

    // APIs
    void run(); // global_place

    void setPlotMode(std::string plot_path);
    void setSolver(const std::string& opt_type);
    void setInitialPlaceMethod(const std::string& init_method);
    void setInitialPlaceOnly();

    // Hyper-Parameter Setting command
    void setTargetOverflow  (float val);
    void setInitLambda      (float val);
    void setInitGammaInv    (float val);
    void setMaxPhiCoef      (float val);
    void setRefHpwl         (float val);
    void setAdamAlpha       (float val);
    void setAdamBeta1       (float val);
    void setAdamBeta2       (float val);
    void setTargetDensity   (float density);

  private:

    // dbDatabase from SkyLine Core
    std::shared_ptr<db::dbDatabase> dbDatabase_;

    // Hyper-Parameters
    std::shared_ptr<HyperParam> param_;

    // 1. Import dbDatabase to SkyPlaceDB
    // 2. Make Sub Tools 
    // (InitialPlacer / Density / WireLength / TargetFunction)
    // 3. Construct TargetFunction f (WireLength & Density) 
    void preamble();

    // Sub-Tools
    std::shared_ptr<SkyPlaceDB>     db_;
    std::shared_ptr<TargetFunction> global_place_problem_;
    std::unique_ptr<InitialPlacer>  initial_placer_;
    std::unique_ptr<SolverBase>     nonlinear_solver_;

    bool localLambdaMode_;
    bool initial_place_only_;

    // Util function
    void freeSubModules();

    // RunTime
    double dbTime_;
    double totalTime_;
};

} // namespace skyplace 

#endif
