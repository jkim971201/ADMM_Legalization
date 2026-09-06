#include <vector>
#include <chrono>
#include <random>    // For mt19937
#include <algorithm> // For sort
#include <cstdio>

#include "object/GPObject.h"
#include "DensityGradient.h"
#include "Util.h"
#include "cuda_linalg/CudaVectorAlgebra.h"

namespace skyplace
{

__device__ inline float getOverlapLength(float cell_min, float cell_max, float bin_min, float bin_max)
{
  return max(0.0f, min(cell_max, bin_max) - max(cell_min, bin_min));
}

__device__ inline OVBIN findBinWithDensitySize(
  const int   numBinX,
  const int   numBinY,
  const float binWidth,
  const float binHeight,
  const float dieLx,
  const float dieLy,
  const float cellCx,
  const float cellCy,
  const float cellDx,
  const float cellDy)
{
  OVBIN ovBins;

  float lx = cellCx - cellDx / 2;
  float ux = cellCx + cellDx / 2;
  
  int minX = floor((lx - dieLx) / binWidth);
  int maxX = ceil ((ux - dieLx) / binWidth);

  minX = max(minX, 0);
  maxX = min(numBinX, maxX);

  ovBins.lxID = minX;
  ovBins.uxID = maxX;

  float ly = cellCy - cellDy / 2;
  float uy = cellCy + cellDy / 2;

  int minY = floor((ly - dieLy) / binHeight);
  int maxY = ceil ((uy - dieLy) / binHeight);

  minY = max(minY, 0);
  maxY = min(numBinY, maxY);

  ovBins.lyID = minY;
  ovBins.uyID = maxY;

  return ovBins;
}

__global__ void computeDensityContributionOfEachCell(
    const int      numMovable,
    const int      numBinX,
    const int      numBinY,
    const float    binWidth,
    const float    binHeight,
    const float    targetDensity,
    const float*   binLx,
    const float*   binLy,
    const float*   binUx,
    const float*   binUy,
    const float    dieLx,
    const float    dieLy,
    const float*   cellCx,
    const float*   cellCy,
    const float*   cellDx,
    const float*   cellDy,
    const float*   cellDensityScale,
    const int*     isFiller,
    const int*     isMacro,
          float*   fillerArea,
          float*   movableArea)
{
  // i := cell_id
  const unsigned int i = blockIdx.x * blockDim.z + threadIdx.z;

  if(i < numMovable)
  {
    bool is_filler = isFiller[i];
    bool is_macro  = isMacro[i];
    float scale = (is_macro == true) ? targetDensity * cellDensityScale[i] : cellDensityScale[i];
    // Macro should be scaled-down with target-density
    // (refer to comments in the OpenROAD RePlAce)

    float cell_lx = cellCx[i] - cellDx[i] / 2.0;
    float cell_ux = cellCx[i] + cellDx[i] / 2.0;
    float cell_ly = cellCy[i] - cellDy[i] / 2.0;
    float cell_uy = cellCy[i] + cellDy[i] / 2.0;

    OVBIN ov_bins
      = findBinWithDensitySize(numBinX, 
                               numBinY,
                               binWidth,
                               binHeight, 
                               dieLx, 
                               dieLy, 
                               cellCx[i],
                               cellCy[i],
                               cellDx[i],
                               cellDy[i]);

    int ov_bins_lx_id = ov_bins.lxID;
    int ov_bins_ly_id = ov_bins.lyID;
    int ov_bins_ux_id = ov_bins.uxID;
    int ov_bins_uy_id = ov_bins.uyID;

    for(int j = ov_bins_lx_id + threadIdx.y; j < ov_bins_ux_id; j += blockDim.y)
    {
      float ov_bin_lx = binLx[j];
      float ov_bin_ux = binUx[j];
      float ov_len_x = getOverlapLength(cell_lx, cell_ux, ov_bin_lx, ov_bin_ux);
      for(int k = ov_bins_ly_id + threadIdx.x; k < ov_bins_uy_id; k+= blockDim.x)
      {
        float ov_bin_ly = binLy[k];
        float ov_bin_uy = binUy[k];
        float ov_len_y = getOverlapLength(cell_ly, cell_uy, ov_bin_ly, ov_bin_uy);

        float ov_area = ov_len_x * ov_len_y * scale;

        int bin_id = j + k * numBinX;
        if(is_filler == true)
          atomicAdd(&(fillerArea[bin_id]), ov_area);
        else
          atomicAdd(&(movableArea[bin_id]), ov_area);
      }
    }
  }
}

__global__ void computeOverflow(
  const int    totalNumBin, 
  const float* fixedArea,
  const float* fillerArea,
  const float* movableArea,
  const float* scaledBinArea,
        float* binDensity,
        float* overflowArea)
{
  const unsigned int binID = blockIdx.x * blockDim.x + threadIdx.x;
  if(binID < totalNumBin)
  {
    float binArea = scaledBinArea[binID];
    binDensity[binID] = (movableArea[binID] + fixedArea[binID] + fillerArea[binID]) / binArea;
    overflowArea[binID] = max(0.0f, movableArea[binID] + fixedArea[binID] - binArea);
  }
}

__global__ void computeMacroOverflow(
  const int    totalNumBin, 
  const float* macroArea,
  const float* scaledBinArea,
        float* macroOverflow)
{
  const unsigned int binID = blockIdx.x * blockDim.x + threadIdx.x;

  if(binID < totalNumBin)
  {
    float binArea = scaledBinArea[binID];
    macroOverflow[binID] = max(0.0f, macroArea[binID] - binArea);
  }
}

__global__ void getDensityGradForEachCell(
  const int    numMovable, 
  const int    numBinX,
  const int    numBinY,
  const float  binWidth,
  const float  binHeight,
  const float* binLx,
  const float* binLy,
  const float* binUx,
  const float* binUy,
  const float  dieLx,
  const float  dieLy,
  const float* cellCx,
  const float* cellCy,
  const float* cellDx,
  const float* cellDy,
  const float* cellDensityScale,
  const float* binLambda,
  const float* electroForceX,
  const float* electroForceY,
        float* densityGradX,
        float* densityGradY,
        float* densityPreconditioner)
{
  // i := cell_id
  const unsigned int i = blockIdx.x * blockDim.z + threadIdx.z;

  if(i < numMovable)
  {
    float cell_lx = cellCx[i] - cellDx[i] / 2.0;
    float cell_ux = cellCx[i] + cellDx[i] / 2.0;
    float cell_ly = cellCy[i] - cellDy[i] / 2.0;
    float cell_uy = cellCy[i] + cellDy[i] / 2.0;

    OVBIN ov_bins
      = findBinWithDensitySize(numBinX, 
                               numBinY,
                               binWidth,
                               binHeight, 
                               dieLx, 
                               dieLy, 
                               cellCx[i],
                               cellCy[i],
                               cellDx[i],
                               cellDy[i]);

    int ov_bins_lx_id = ov_bins.lxID;
    int ov_bins_ly_id = ov_bins.lyID;
    int ov_bins_ux_id = ov_bins.uxID;
    int ov_bins_uy_id = ov_bins.uyID;

    extern __shared__ float shared[]; 
    // size will be sizeof(float) * blockDim.z * 3

    float* shared_x = shared;                // this is for x-direction force
    float* shared_y = shared   + blockDim.z; // this is for y-direction force
    float* shared_p = shared_y + blockDim.z; // this is for density preconditioner

    if(threadIdx.x == 0 && threadIdx.y == 0) // Master Thread for this cell
    {
      // if threadIdx.z are same (in each block), the kernel is updating the same cell.
      shared_x[threadIdx.z] = 0.0f;
      shared_y[threadIdx.z] = 0.0f;
      shared_p[threadIdx.z] = 0.0f;
    }
    __syncthreads();

    float force_x_this_thread = 0.0f;
    float force_y_this_thread = 0.0f;
    float precond_this_thread = 0.0f;

    for(int j = ov_bins_lx_id + threadIdx.y; j < ov_bins_ux_id; j += blockDim.y)
    {
      float ov_bin_lx = binLx[j];
      float ov_bin_ux = binUx[j];
      float ov_len_x = getOverlapLength(cell_lx, cell_ux, ov_bin_lx, ov_bin_ux);
      for(int k = ov_bins_ly_id + threadIdx.x; k < ov_bins_uy_id; k+= blockDim.x)
      {
        float ov_bin_ly = binLy[k];
        float ov_bin_uy = binUy[k];
        float ov_len_y = getOverlapLength(cell_ly, cell_uy, ov_bin_ly, ov_bin_uy);

        float ov_area = ov_len_x * ov_len_y * cellDensityScale[i];

        int bin_id = j + k * numBinX;
        float lambda = binLambda[bin_id]; // local lambda 

        // GPU version of Poisson Solver is
        // multiplied by sqrt(2) * sqrt(2) 
        // Therefore, electroForce is multiplied by 0.5
        force_x_this_thread += lambda * ov_area * electroForceX[bin_id] * 0.5;
        force_y_this_thread += lambda * ov_area * electroForceY[bin_id] * 0.5;
        precond_this_thread += lambda * ov_area;
      }
    }

    atomicAdd(&shared_x[threadIdx.z], force_x_this_thread);
    atomicAdd(&shared_y[threadIdx.z], force_y_this_thread);
    atomicAdd(&shared_p[threadIdx.z], precond_this_thread);
    __syncthreads();

    if(threadIdx.x == 0 && threadIdx.y == 0) // Master Thread for this cell
    {
      densityGradX[i] = shared_x[threadIdx.z];
      densityGradY[i] = shared_y[threadIdx.z];
      densityPreconditioner[i] = shared_p[threadIdx.z];
    }
  }
}

__global__ void updateLocalLambda(
  const int    totalNumBin, 
  const float  targetDensity,
  const float* binDensity,
        float* binLambda)
{
  const unsigned int binID = blockIdx.x * blockDim.x + threadIdx.x;

  if(binID < totalNumBin)
  {
    float density = binDensity[binID];

    if(density > targetDensity)
      binLambda[binID] = __powf(1 + __log10f(density), 0.5);
    else
      binLambda[binID] = 1.0;
  }
}

void
DensityGradient::computeGrad(
       float* densityGradX,
       float* densityGradY,
 const float* cellCx,
 const float* cellCy)
{
  int numThreadBin = 256;
  int numBlockBin  = (numBinX_ * numBinY_ - 1 + numThreadBin) / numThreadBin;

  auto density_start = getChronoNow();

  const int  z_size_per_block = 64;
  const int  cell_kernel_num_block = (numMovable_ + z_size_per_block - 1) / z_size_per_block;
  const dim3 cell_kernel_block_size(2, 2, z_size_per_block);
  size_t shared_mem = sizeof(float) * z_size_per_block * 3;
  // shared_mem := {X-direction force, Y-direction force, preconditioner}

  cudaDeviceSynchronize();

  // Step 1. Initialize movableArea as zero
  d_movableArea_.fillZero();
  d_fillerArea_.fillZero();

  // Step 2. Compute OverlapArea
  computeDensityContributionOfEachCell<<<cell_kernel_num_block, cell_kernel_block_size>>>(
    numMovable_, 
    numBinX_, 
    numBinY_, 
    binWidth_, 
    binHeight_,
    targetDensity_,
    d_ptr_binLx_,
    d_ptr_binLy_,
    d_ptr_binUx_,
    d_ptr_binUy_,
    dieLx_, 
    dieLy_, 
    cellCx,
    cellCy,
    d_ptr_cellDensityWidth_,
    d_ptr_cellDensityHeight_,
    d_ptr_cellDensityScale_,
    d_ptr_isFiller_,
    d_ptr_isMacro_,
    d_ptr_fillerArea_,
    d_ptr_movableArea_);

  // Step 3. Update Overflow
  computeOverflow<<<numBlockBin, numThreadBin>>>(
    numBinX_ * numBinY_,
    d_ptr_fixedArea_,
    d_ptr_fillerArea_,
    d_ptr_movableArea_,
    d_ptr_scaledBinArea_,
    d_ptr_binDensity_,
    d_ptr_overflowArea_);

  overflow_ = computeVectorSum(d_overflowArea_) / sumMovableArea_;

  cudaDeviceSynchronize();
  auto overflow_finish = getChronoNow();
  binDenUpTime_ += evalTime(density_start);

  // Step 4. Solve Poisson Equation
  poissonSolver_->solvePoisson(d_ptr_binDensity_, 
                               d_ptr_binPotential_, 
                               d_ptr_electroForceX_, 
                               d_ptr_electroForceY_);

  cudaDeviceSynchronize();
  auto poisson_finish = getChronoNow();
  poissonTime_ += evalTime(overflow_finish);

  // Step 5. Compute Density Gradient based on Electric Force XY
  getDensityGradForEachCell<<<cell_kernel_num_block, cell_kernel_block_size, shared_mem>>>(
    numMovable_, 
    numBinX_,
    numBinY_,
    binWidth_, 
    binHeight_,
    d_ptr_binLx_,
    d_ptr_binLy_,
    d_ptr_binUx_,
    d_ptr_binUy_,
    dieLx_, 
    dieLy_,
    cellCx,
    cellCy,
    d_ptr_cellDensityWidth_,
    d_ptr_cellDensityHeight_,
    d_ptr_cellDensityScale_,
    d_ptr_binLambda_,
    d_ptr_electroForceX_,
    d_ptr_electroForceY_,
    densityGradX,
    densityGradY,
    d_ptr_densityPreconditioner_);

  cudaDeviceSynchronize();
  densityTime_ += evalTime(density_start);
}

void
DensityGradient::initForCUDAKernel()
{
  const int numBin2 = numBinX_ * numBinY_;
  std::vector<float> h_binLambda(numBin2);
  std::vector<float> h_fixedArea(numBin2);
  std::vector<float> h_scaledBinArea(numBin2);
  std::vector<float> h_binLx(numBinX_);
  std::vector<float> h_binUx(numBinX_);
  std::vector<float> h_binLy(numBinY_);
  std::vector<float> h_binUy(numBinY_);

  d_isFiller_.resize(numMovable_);
  d_isMacro_.resize(numMovable_);

  d_fixedArea_.resize(numBin2);
  d_macroArea_.resize(numBin2);
  d_scaledBinArea_.resize(numBin2);

  d_binLx_.resize(numBinX_);
  d_binUx_.resize(numBinX_);
  d_binLy_.resize(numBinY_);
  d_binUy_.resize(numBinY_);

  d_movableArea_.resize(numBin2);
  d_fillerArea_.resize(numBin2);
  d_overflowArea_.resize(numBin2);
  d_macroOverflowArea_.resize(numBin2);

  d_binDensity_.resize(numBin2);
  d_binLambda_.resize(numBin2);
  d_binPotential_.resize(numBin2);
  d_binPenalty_.resize(numBin2);

  d_densityPreconditioner_.resize(numMovable_);

  d_electroForceX_.resize(numBin2);
  d_electroForceY_.resize(numBin2);

  d_cellDensityWidth_.resize(numMovable_);
  d_cellDensityHeight_.resize(numMovable_);
  d_cellDensityScale_.resize(numMovable_);

  d_ptr_isFiller_       = d_isFiller_.data();
  d_ptr_isMacro_        = d_isMacro_.data();

  d_ptr_fixedArea_      = d_fixedArea_.data();
  d_ptr_macroArea_      = d_macroArea_.data();
  d_ptr_scaledBinArea_  = d_scaledBinArea_.data();

  d_ptr_binLx_          = d_binLx_.data();
  d_ptr_binUx_          = d_binUx_.data();
  d_ptr_binLy_          = d_binLy_.data();
  d_ptr_binUy_          = d_binUy_.data();

  d_ptr_movableArea_       = d_movableArea_.data();
  d_ptr_fillerArea_        = d_fillerArea_.data();
  d_ptr_overflowArea_      = d_overflowArea_.data();
  d_ptr_macroOverflowArea_ = d_macroOverflowArea_.data();

  d_ptr_binDensity_     = d_binDensity_.data();
  d_ptr_binLambda_      = d_binLambda_.data();
  d_ptr_binPotential_   = d_binPotential_.data();
  d_ptr_binPenalty_     = d_binPenalty_.data();

  d_ptr_densityPreconditioner_ = d_densityPreconditioner_.data();

  d_ptr_electroForceX_  = d_electroForceX_.data();
  d_ptr_electroForceY_  = d_electroForceY_.data();

  d_ptr_cellDensityWidth_  = d_cellDensityWidth_.data();
  d_ptr_cellDensityHeight_ = d_cellDensityHeight_.data();
  d_ptr_cellDensityScale_  = d_cellDensityScale_.data();

  const auto& db_bins = db_->bins();
  for(int bin_id = 0; bin_id < numBin2; bin_id++)
  {
    auto bin = db_bins[bin_id];
    h_fixedArea[bin_id] = bin->fixedArea();
    h_scaledBinArea[bin_id] = bin->area() * bin->targetDensity();
    h_binLambda[bin_id] = bin->lambda();
  }

  for(int col = 0; col < numBinX_; col++)
  {
    int bin_id = col;
    auto bin = db_bins[bin_id];
    // Assume every bin[x][*] has same width 
    h_binLx[col] = bin->lx();
    h_binUx[col] = bin->ux();
  }

  for(int row = 0; row < numBinY_; row++)
  {
    int bin_id = row * numBinX_;
    auto bin = db_bins[bin_id];
    // Assume every bin[*][y] has same height
    h_binLy[row] = bin->ly();
    h_binUy[row] = bin->uy();
  }

  std::vector<int>   h_isFiller(numMovable_);
  std::vector<int>   h_isMacro(numMovable_);
  std::vector<float> h_cellDensityWidth(numMovable_);
  std::vector<float> h_cellDensityHeight(numMovable_);
  std::vector<float> h_cellDensityScale(numMovable_);

  int i = 0;
  for(auto& cell : db_->movableCells())
  {
    h_isFiller[i]          = cell->isFiller();
    h_isMacro[i]           = cell->isMacro();
    h_cellDensityWidth[i]  = cell->dDx();
    h_cellDensityHeight[i] = cell->dDy();
    h_cellDensityScale[i]  = cell->densityScale();
    i++;
  }

  // Host -> Device
  d_isFiller_          = h_isFiller;
  d_isMacro_           = h_isMacro;
  d_cellDensityWidth_  = h_cellDensityWidth;
  d_cellDensityHeight_ = h_cellDensityHeight;
  d_cellDensityScale_  = h_cellDensityScale;
  d_fixedArea_         = h_fixedArea;
  d_scaledBinArea_     = h_scaledBinArea;
  d_binLx_             = h_binLx;
  d_binLy_             = h_binLy;
  d_binUx_             = h_binUx;
  d_binUy_             = h_binUy;
  d_binLambda_         = h_binLambda;
}

DensityGradient::DensityGradient()
  : numBinX_         (0),
    numBinY_         (0),
    numMovable_      (0),
    sumMovableArea_  (0),
    overflow_        (0),
    dieLx_           (0),
    dieLy_           (0),
    dieUx_           (0),
    dieUy_           (0),
    targetDensity_   (0),
    binWidth_        (0),
    binHeight_       (0),
    localLambdaMode_ (false),

    densityTime_     (0.0),
    poissonTime_     (0.0),
    binDenUpTime_    (0.0)
{}

DensityGradient::DensityGradient(std::shared_ptr<SkyPlaceDB> db)
  : DensityGradient()
{
  db_ = db;

  numBinX_        = db_->numBinX();
  numBinY_        = db_->numBinY();
  numMovable_     = db_->numMovable();
  sumMovableArea_ = db_->sumMovableArea();

  dieLx_ = db_->die()->lx();
  dieLy_ = db_->die()->ly();
  dieUx_ = db_->die()->ux();
  dieUy_ = db_->die()->uy();

  binWidth_  = db_->binX();
  binHeight_ = db_->binY();

  targetDensity_ = db_->targetDensity();
  poissonSolver_ = std::make_unique<PoissonSolver>(numBinX_, numBinY_);

  initForCUDAKernel();
}

} // namespace skyplace
