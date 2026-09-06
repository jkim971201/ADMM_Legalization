#include "cuda_linalg/CudaVectorAlgebra.h"

#include "Util.h"
#include "HyperParam.h"
#include "SkyPlaceDB.h"
#include "object/GPObject.h"
#include "wirelength_gradient/WireLengthGradient.h"
#include "density_gradient/DensityGradient.h"

#include "TargetFunction.h"
#include "ProblemInstanceKernel.hpp"

namespace skyplace
{

TargetFunction::TargetFunction(
  std::shared_ptr<SkyPlaceDB> db, 
  std::shared_ptr<HyperParam> param)
  : db_(db), param_(param)
{
  // Initialize Hyperparameters
  lambda_     = param->initLambda;
  minPrecond_ = param->minPrecond;
  gammaInv_   = 2.0 * param->initGammaInvCoef / (db_->binX() + db_->binY());

  x_min_ = db_->die()->lx();
  y_min_ = db_->die()->ly();
  x_max_ = db_->die()->ux();
  y_max_ = db_->die()->uy();

  // Make SubModules
  wireLength_ = std::make_shared<WireLengthGradient>(db_);
  wireLength_->setGammaInv(gammaInv_);

  density_ = std::make_shared<DensityGradient>(db_);

  if(param_->plot_mode == true)
    painter_ = std::make_unique<Painter>(db_);

  prepareData();
}

void
TargetFunction::prepareData()
{
  auto& movable_cells = db_->movableCells();

  num_movable_ = movable_cells.size();
  num_var_ = 2 * num_movable_;

  // Host-side
  h_cell_pos_.resize(num_var_);

  std::vector<float> cell2num_pin;
  std::vector<float> cell_real_width;
  std::vector<float> cell_real_height;

	cell2num_pin.reserve(db_->numPin());
  cell_real_width.reserve(num_movable_);
  cell_real_height.reserve(num_movable_);

  // Device-side
  d_wl_grad_.resize(num_var_, 0.0);
  d_density_grad_.resize(num_var_, 0.0);

  d_cell2num_pin_.resize(num_movable_, 0.0);
  d_cell_real_width_.resize(num_movable_);
  d_cell_real_height_.resize(num_movable_);

  int cell_id = 0;
  for(auto& cell : movable_cells)
  {
    cell2num_pin.push_back(static_cast<float>(cell->pins().size()));
    cell_real_width.push_back(cell->dx());
    cell_real_height.push_back(cell->dy());
    h_cell_pos_[cell_id] = cell->cx();
    h_cell_pos_[cell_id + num_movable_] = cell->cy();
    cell_id++;
  }

  // Host to Device
  d_cell2num_pin_     = cell2num_pin;
  d_cell_real_width_  = cell_real_width;
  d_cell_real_height_ = cell_real_height;
}

void
TargetFunction::updatePointAndGetGrad(const CudaVector<float>& var, CudaVector<float>& grad)
{
  const float* x_var = var.data();
  const float* y_var = var.data() + num_movable_;

  float* wl_grad_x = d_wl_grad_.data();
  float* wl_grad_y = d_wl_grad_.data() + num_movable_;

  float* den_grad_x = d_density_grad_.data();
  float* den_grad_y = d_density_grad_.data() + num_movable_;

  float* grad_x = grad.data();
  float* grad_y = grad.data() + num_movable_;

  // For WireLength Grad
  wireLength_->updatePinCoordinates(x_var, y_var);

  wireLength_->computeGrad(wl_grad_x, wl_grad_y);

  // For Density Grad
  density_->computeGrad(den_grad_x, den_grad_y, x_var, y_var);

  int num_thread = 64;
  int num_block_cell = (num_movable_ - 1 + num_thread) / num_thread;

  bool precondition = true;
  if(precondition == true)
  {
    jacobianPrecondition<<<num_block_cell, num_thread>>>(
      num_movable_, 
      lambda_,
      minPrecond_,
      d_cell2num_pin_.data(),
      density_->getDevicePreconditioner(),
      wl_grad_x,
      wl_grad_y,
      den_grad_x,
      den_grad_y,
      grad_x,
      grad_y);
  }
  else
  {
    noPrecondition<<<num_block_cell, num_thread>>>(
      num_movable_, 
      lambda_,
      wl_grad_x,
      wl_grad_y,
      den_grad_x,
      den_grad_y,
      grad_x,
      grad_y);
  }
}

void
TargetFunction::clipToChipBoundary(CudaVector<float>& cell_pos)
{
  int num_thread = 64;
  int num_block_cell = (num_movable_ - 1 + num_thread) / num_thread;

  clipToChipBoundaryKernel<<<num_block_cell, num_thread>>>(
    num_movable_,
    x_min_,
    y_min_,
    x_max_,
    y_max_,
    d_cell_real_width_.data(),
    d_cell_real_height_.data(),
    cell_pos.data(),
    cell_pos.data() + num_movable_);
}

void
TargetFunction::importFromDb()
{
  auto& db_movable_cells = db_->movableCells();
  int num_movable = db_movable_cells.size();
  for(int cell_id = 0; cell_id < num_movable; cell_id++)
  {
    GPCell* cell = db_movable_cells[cell_id];
    h_cell_pos_[cell_id] = cell->cx();
    h_cell_pos_[cell_id + num_movable] = cell->cy();
  }
}

void
TargetFunction::exportToSolver(CudaVector<float>& var_from_solver)
{
  // Host to Device
  var_from_solver = h_cell_pos_;
}

void
TargetFunction::getInitialGrad(
  const CudaVector<float>& initial_var,
        CudaVector<float>& initial_grad)
{
  const float* x_var = initial_var.data();
  const float* y_var = initial_var.data() + num_movable_;

  float* wl_grad_x = d_wl_grad_.data();
  float* wl_grad_y = d_wl_grad_.data() + num_movable_;

  float* den_grad_x = d_density_grad_.data();
  float* den_grad_y = d_density_grad_.data() + num_movable_;

  float* grad_x = initial_grad.data();
  float* grad_y = initial_grad.data() + num_movable_;

  // We have to first compute Density Gradient
  // because we do not yet know overflow that is required to update gamma
  density_->computeGrad(den_grad_x, den_grad_y, x_var, y_var);

  overflow_ = density_->overflow();
  updateGammaInv(overflow_);

  wireLength_->updatePinCoordinates(x_var, y_var);
  wireLength_->computeGrad(wl_grad_x, wl_grad_y);

  int num_thread = 64;
  int num_block_cell = (num_movable_ - 1 + num_thread) / num_thread;

  jacobianPrecondition<<<num_block_cell, num_thread>>>(
    num_movable_, 
    0.0,
    minPrecond_,
    d_cell2num_pin_.data(),
    density_->getDevicePreconditioner(),
    wl_grad_x,
    wl_grad_y,
    den_grad_x,
    den_grad_y,
    grad_x,
    grad_y);

  float wlGradSum = computeNorm1(d_wl_grad_);
  float densityGradSum = computeNorm1(d_density_grad_);

  if(std::isnan(wlGradSum) || std::isnan(densityGradSum) ||
     std::isinf(wlGradSum) || std::isinf(densityGradSum))
  {
    printf("DensityGradSum   : %f\n", densityGradSum);
    printf("WireLengthGradSum: %f\n", wlGradSum);
    exit(0);
  }

  lambda_ *= wlGradSum / densityGradSum;
  // Set the initial value of lambda
  // same as RePlAce Paper TCAD 2019
}

void
TargetFunction::updateGammaInv(float overflow)
{
  if(overflow > 1.0)
    gammaInv_ = 0.1;
  else if(overflow < 0.1)
    gammaInv_ = 10.0;
  else
    gammaInv_ = 1.0 / std::pow(10.0, 20 * (overflow - 0.1) / 9.0 - 1.0);

  gammaInv_ *= param_->initGammaInvCoef;
  // 10, 20, 9 are also hyper-parameters...

  wireLength_->setGammaInv(gammaInv_);
}

void
TargetFunction::updateLambda(float prevHpwl, float curHpwl)
{
  float p = (curHpwl - prevHpwl) / param_->referenceHpwl;
  float coef = 1;

  // See ePlace paper for this udpate scheme
  if(p < 0) 
    coef = param_->maxPhiCoef;
  else
  {
    coef 
      = std::max(param_->minPhiCoef, 
                 param_->maxPhiCoef * std::pow(param_->maxPhiCoef, -p));
  }
  lambda_ *= coef;
}

void
TargetFunction::updateParameters()
{
  prevHpwl_ = hpwl_;
  hpwl_     = wireLength_->computeHPWL();
  overflow_ = density_->overflow();

  updateGammaInv(overflow_);
  updateLambda(prevHpwl_, hpwl_);
}

bool
TargetFunction::checkConvergence() const
{
  bool convergence = false;
  if(overflow_ <= param_->targetOverflow)
  {
    // convergence is detected but keep going
    // if hpwl is descreasing
    if(hpwl_ < prevHpwl_)
      convergence = false; 
    else
      convergence = true;
  }
  return convergence;
}

int
TargetFunction::getNumVariable() const
{
  return num_var_;
}

void 
TargetFunction::solveBgnCbk() 
{
  printf("Initial HPWL     : %-6.3f\n", hpwl_ / 1e6);
  printf("Initial Overflow : %-3.3f\n", overflow_);
  printMetricRow();

  // Prepare Painting
  if(param_->plot_mode == true)
    painter_->preparePlot(param_->plot_path);
}

void
TargetFunction::solveEndCbk(int iter, double runtime, const CudaVector<float>& var)
{
  printProgress(iter, runtime); // Last progress message
  
  printf("----------------------------------------------\n");
  printf(" Design   | %-10s  \n", db_->designName().c_str());
  printf(" HPWL     | %-10.5f\n", hpwl_);
  printf(" Iter     | %-10d  \n", iter);
  printf(" Overflow | %-10.5f\n", overflow_);
  printf("----------------------------------------------\n");

  exportToDb(var);

  if(param_->plot_mode == true)
    painter_->saveImage(iter, hpwl_, overflow_); 
}

void
TargetFunction::iterBgnCbk(int iter)
{
  // Empty Callback
}

void
TargetFunction::iterEndCbk(int iter, double runtime, const CudaVector<float>& var)
{
  updateParameters(); 
  if(iter == 0 || iter % param_->log_freq == 0)
    printProgress(iter, runtime);
  // This should go to iterBgnCbk?

  if(param_->plot_mode == true && (iter % param_->plot_freq == 0))
  {
    exportToDb(var);
    painter_->saveImage(iter, hpwl_, overflow_);
  }
}

void
TargetFunction::printMetricRow() const
{
  printf("\n");
  printf("----------------------------------------------\n");
  printf("|  Iter |     HPWL     | Overflow | Time (s) |\n");
  printf("----------------------------------------------\n");
}

void
TargetFunction::printProgress(int iter, double runtime) const
{
  printf("| %5d | %12.1f |  %-4.4f  | %8.2f |\n",
          iter, hpwl_, overflow_, runtime);
}

void
TargetFunction::exportToDb(const CudaVector<float>& cell_pos)
{
  thrust::copy(cell_pos.begin(), cell_pos.end(), h_cell_pos_.begin());

  auto& db_movable_cells = db_->movableCells();
  int num_movable = db_movable_cells.size();
  for(int cell_id = 0; cell_id < num_movable; cell_id++)
  {
    GPCell* cell = db_movable_cells[cell_id];
    cell->setCenterLocation(h_cell_pos_[cell_id], h_cell_pos_[cell_id + num_movable]);
  }
}

void
TargetFunction::diverge()
{
  printf("Divergence Detected!\n");
  printf("Terminate Placer...\n");
  exit(0);
}

} // namespace skyplace
