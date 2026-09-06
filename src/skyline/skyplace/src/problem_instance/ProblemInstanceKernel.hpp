#ifndef PROBLEM_INSTANCE_KERNEL_H
#define PROBLEM_INSTANCE_KERNEL_H

namespace skyplace
{

__device__ inline float getXInsideChip(
  const float cell_cx, 
  const float cell_width, 
  const float x_min,
  const float x_max)
{
  float new_cx = cell_cx;
  if(cell_cx - cell_width / 2 < x_min)
    new_cx = x_min + cell_width / 2;
  if(cell_cx + cell_width / 2 > x_max)
    new_cx = x_max - cell_width / 2;
  return new_cx;
}

__device__ inline float getYInsideChip(
  const float cell_cy, 
  const float cell_height, 
  const float y_min,  
  const float y_max)
{
  float new_cy = cell_cy;
  if(cell_cy - cell_height / 2 < y_min)
    new_cy = y_min + cell_height / 2;
  if(cell_cy + cell_height / 2 > y_max)
    new_cy = y_max - cell_height / 2;
  return new_cy;
}

__global__ void clipToChipBoundaryKernel(
  const int    num_cell,
  const float  x_min,
  const float  y_min,
  const float  x_max,
  const float  y_max,
  const float* cell_width,
  const float* cell_height,
        float* cell_cx,
        float* cell_cy)
{
  const int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(cell_id < num_cell)
  {
    cell_cx[cell_id] = getXInsideChip(cell_cx[cell_id], cell_width[cell_id], x_min, x_max);
    cell_cy[cell_id] = getYInsideChip(cell_cy[cell_id], cell_height[cell_id], y_min, y_max);
  }
}

__global__ void noPrecondition(
  const int    numCell,
  const float  lambda,
  const float* wlGradX,
  const float* wlGradY,
  const float* densityGradX,
  const float* densityGradY,
        float* totalGradX,
        float* totalGradY)
{
  // i := cellID
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < numCell)
  {
    totalGradX[i] = (wlGradX[i] + lambda * densityGradX[i]);
    totalGradY[i] = (wlGradY[i] + lambda * densityGradY[i]);
  }
}

__global__ void jacobianPrecondition(
  const int    numCell,
  const float  lambda,
  const float  minPrecond,
  const float* numPinForEachCell,
  const float* densityPreconditioner,
  const float* wlGradX,
  const float* wlGradY,
  const float* densityGradX,
  const float* densityGradY,
        float* totalGradX,
        float* totalGradY)
{
  // i := cellID
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < numCell)
  {
    int   wirelengthPrecond = numPinForEachCell[i];
    float densityPrecond    = densityPreconditioner[i];
    float precond = wirelengthPrecond + lambda * densityPrecond;
    precond = max(precond, minPrecond);

	  totalGradX[i] = (wlGradX[i] + lambda * densityGradX[i]) / precond;
	  totalGradY[i] = (wlGradY[i] + lambda * densityGradY[i]) / precond;
  }
}

}

#endif
