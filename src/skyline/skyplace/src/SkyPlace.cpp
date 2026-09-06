#include <cstdio>
#include <chrono>

#include "skyplace/SkyPlace.h"
#include "SkyPlaceDB.h"

#include "db/dbDatabase.h"
#include "db/dbTech.h"
#include "db/dbDesign.h"

#include "object/GPObject.h"

#include "Util.h"
#include "HyperParam.h"
#include "problem_instance/TargetFunction.h"
#include "initial_place/InitialPlacer.h"
#include "nonlinear_solver/AdamSolver.h"
#include "nonlinear_solver/NesterovSolver.h"

namespace skyplace 
{

SkyPlace::SkyPlace() {}
SkyPlace::SkyPlace(std::shared_ptr<dbDatabase> db)
  : dbDatabase_         (db),
    localLambdaMode_    (false),
    initial_place_only_ (false)
{
  db_    = std::make_shared<SkyPlaceDB>();
  param_ = std::make_shared<HyperParam>();
}

SkyPlace::~SkyPlace() {}

void
SkyPlace::setPlotMode(std::string plot_path)
{
  param_->plot_mode = true;
  param_->plot_path = plot_path;
}

void
SkyPlace::setInitialPlaceOnly()
{
  initial_place_only_ = true;
}

void 
SkyPlace::setTargetOverflow(float val) 
{ 
  param_->targetOverflow = val;      
}

void 
SkyPlace::setInitLambda(float val) 
{ 
  param_->initLambda = val;
}

void 
SkyPlace::setMaxPhiCoef(float val) 
{ 
  param_->maxPhiCoef = val;
}

void 
SkyPlace::setRefHpwl(float val) 
{ 
  param_->referenceHpwl = val;
}

void 
SkyPlace::setInitGammaInv(float val) 
{ 
  param_->initGammaInvCoef = val;
}

void 
SkyPlace::setAdamAlpha(float val) 
{ 
  param_->adam_alpha = val;           
}

void 
SkyPlace::setAdamBeta1(float val) 
{ 
  param_->adam_beta1 = val;           
}

void 
SkyPlace::setAdamBeta2(float val) 
{ 
  param_->adam_beta2 = val;           
}

// Though targetDensity is also a hyper-parameter,
// this has to be initialized before SkyPlaceDB is initialized.
void 
SkyPlace::setTargetDensity(float density)  
{ 
  float target_density_percentage         = density * 100.0;
  float target_density_default            = 1.0;
  float target_density_default_percentage = target_density_default * 100.0;

  if(target_density_percentage > 100.00 || target_density_percentage <= 0.0)
  {
    printf("Invalid target density (%4.2f%%). %4.2f%% will be used.\n", 
        target_density_percentage, target_density_default_percentage);
    db_->setTargetDensity(target_density_default);
  }
  else
    db_->setTargetDensity(density);
}

void
SkyPlace::preamble()
{
  /* SkyPlaceDB Initialization */
  auto preamble_start = getChronoNow();
  db_->init(dbDatabase_);
  // Runtime messages will be printed inside init() function.
  // This is to print runtime of each subprocedures.

  /* InitialPlacer Initialization */
  auto ip_init_start = getChronoNow();
  initial_placer_ = std::make_unique<InitialPlacer>(db_);
  const double ip_init_time = evalTime(ip_init_start);
  printf("Initialize InitialPlacer      (takes %5.2f s)\n", ip_init_time);

  /* Problem Instance Initialization */
  auto problem_init_start = getChronoNow();
  global_place_problem_ = std::make_shared<TargetFunction>(db_, param_);
  const double problem_init_time = evalTime(problem_init_start);
  printf("Initialize GlobalPlaceProblem (takes %5.2f s)\n", problem_init_time);

  /* Solver Initialization */
  auto solver_init_start = getChronoNow();

  if(param_->solver_type == ADAM) 
    nonlinear_solver_ = std::make_unique<AdamSolver>(param_, global_place_problem_);
  else
    nonlinear_solver_ = std::make_unique<NesterovSolver>(param_, global_place_problem_);
                           
  const double solver_init_time = evalTime(solver_init_start);
  printf("Initialize NonlinearSolver    (takes %5.2f s)\n", solver_init_time);

  dbTime_ = evalTime(preamble_start);
  printf("--> Total Initialization takes %5.2f s\n", dbTime_);

  db_->printInfo();
}

void
SkyPlace::setInitialPlaceMethod(const std::string& init_method)
{
  if(init_method == "random")
    param_->init_method = InitialPlaceMethod::RANDOM_START;
  else if(init_method == "sdp")
    param_->init_method = InitialPlaceMethod::SDP_RELAXATION;
  else
  {
    printf("[WARN] Undefined Initialization method %s\n", init_method.c_str());
    printf("[WARN] Argument is ignored and random initialization will be used...\n");
  }
}

void
SkyPlace::setSolver(const std::string& solver_type_str)
{
  if(solver_type_str == "Adam" || solver_type_str == "adam")
    param_->solver_type = SolverType::ADAM;
  else if(solver_type_str == "Nesterov" || solver_type_str == "nesterov")
    param_->solver_type = SolverType::NESTEROV;
  else
  {
    printf("Undefined SolverType %s\n", solver_type_str.c_str());
    printf("Argument is ignored and Nesterov will be used...\n");
  }
}

void
SkyPlace::run()
{
  auto run_start = std::chrono::high_resolution_clock::now();

  printf("\n");
  printf("Run global_place...\n");

  // Prepare for Global Placement
  // (Import dbDatabase to SkyPlaceDB + Make sub-tools...)
  preamble();

  /* Step #1: Initial Placement */
  initial_placer_->doInitialPlace(param_->init_method);

  db_->exportDB(dbDatabase_);            // Deliver new coorindates to dbDatabase

  if(initial_place_only_ == false)
  {
    // Synchronize TargetFunction and SkyPlaceDB
    global_place_problem_->importFromDb();

    /* Step #2: Main Placement Iteration */
    nonlinear_solver_->solve();

     // TODO : Fix mismatch between hpwl_ of Optimizer class and SkyPlaceDB
    db_->updateHpwl();
    db_->exportDB(dbDatabase_);  // Deliver new coorindates to dbDatabase

    totalTime_ = evalTime(run_start);
  }

  freeSubModules();
}

void
SkyPlace::freeSubModules()
{
  // Calling a CUDA API in the constructor or destructor
  // of static variable is undefined behavior.
  // This is because constructor/destructor of global static variable 
  // is called before/after main function.
  // Also, it is unsafe to call CUDA API in the destructor.
  // Therefore, we will free device memory by this way.
  // NOTE : This will cause memory leak in the PoissonSolver.
  // (since it does not use thrust, its memory will not be released automatically)
  // TODO : Fix this memory leak.
  nonlinear_solver_.reset();
}

} // namespace SkyPlace
