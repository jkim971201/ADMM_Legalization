#include <cstdio>
#include <cmath>
#include <chrono>
#include <climits> // for MAX FLOAT

#include "cuda_linalg/CudaVectorAlgebra.h"
#include "object/GPObject.h"
#include "WireLengthGradient.h"
#include "Util.h"

namespace skyplace
{

__global__ void computeMinMax(
  const int    numNet,
  const int*   netStart,
  const float* pinX, 
  const float* pinY, 
        float* maxPinXArr, 
        float* minPinXArr,
        float* maxPinYArr, 
        float* minPinYArr)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < numNet)
  {
    float maxPinX = 0;
    float maxPinY = 0;
    float minPinX = k_float_max;
    float minPinY = k_float_max;

    for(int j = netStart[i]; j < netStart[i+1]; j++)
    {
      maxPinX = max(pinX[j], maxPinX);
      minPinX = min(pinX[j], minPinX);
      maxPinY = max(pinY[j], maxPinY);
      minPinY = min(pinY[j], minPinY);
    }

    maxPinXArr[i] = maxPinX;
    minPinXArr[i] = minPinX;
    maxPinYArr[i] = maxPinY;
    minPinYArr[i] = minPinY;
  }
}

__global__ void computeABCKernel(
  const int    numPin,
  const float  gammaInv, 
  const int*   pin2Net,
  const float* pinX, 
  const float* pinY, 
  const float* maxPinX,
  const float* maxPinY,
  const float* minPinX,
  const float* minPinY,
        float* apX, 
        float* apY, 
        float* amX, 
        float* amY,
        float* bpX, 
        float* bpY, 
        float* bmX, 
        float* bmY,
        float* cpX, 
        float* cpY, 
        float* cmX, 
        float* cmY)
{
  const unsigned int pin_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(pin_id < numPin)
  {
    int net_id = pin2Net[pin_id];

    float pin_x = pinX[pin_id];
    float pin_y = pinY[pin_id];
    float max_pin_x = maxPinX[net_id];
    float max_pin_y = maxPinY[net_id];
    float min_pin_x = minPinX[net_id];
    float min_pin_y = minPinY[net_id];

    float ap_x = exp(+(pin_x - max_pin_x) * gammaInv);  
    float ap_y = exp(+(pin_y - max_pin_y) * gammaInv);  
    float am_x = exp(-(pin_x - min_pin_x) * gammaInv);  
    float am_y = exp(-(pin_y - min_pin_y) * gammaInv);  

    apX[pin_id] = ap_x;
    apY[pin_id] = ap_y;
    amX[pin_id] = am_x;
    amY[pin_id] = am_y;

    atomicAdd(&bpX[net_id], ap_x);
    atomicAdd(&bpY[net_id], ap_y);
    atomicAdd(&bmX[net_id], am_x);
    atomicAdd(&bmY[net_id], am_y);

    atomicAdd(&cpX[net_id], pin_x * ap_x);
    atomicAdd(&cpY[net_id], pin_y * ap_y);
    atomicAdd(&cmX[net_id], pin_x * am_x);
    atomicAdd(&cmY[net_id], pin_y * am_y);
  }
}

__global__ void computePinGrad(
  const int    numPin,
  const float  gInv,
  const int*   pin2Net,
  const float* pinX,
  const float* pinY,
  const float* apX, 
  const float* apY,
  const float* amX, 
  const float* amY,
  const float* bpX, 
  const float* bpY,
  const float* bmX, 
  const float* bmY,
  const float* cpX, 
  const float* cpY,
  const float* cmX, 
  const float* cmY,
  const float* netWeight,
        float* pinGradX,  
        float* pinGradY)
{
  // i := pin_id
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < numPin)
  {
    int nID = pin2Net[i]; // netID
    pinGradX[i] 
      = ((1 - pinX[i] * gInv) * bmX[nID] + gInv * cmX[nID]) * amX[i] / bmX[nID] / bmX[nID]
      - ((1 + pinX[i] * gInv) * bpX[nID] - gInv * cpX[nID]) * apX[i] / bpX[nID] / bpX[nID];

    pinGradY[i] 
      = ((1 - pinY[i] * gInv) * bmY[nID] + gInv * cmY[nID]) * amY[i] / bmY[nID] / bmY[nID]
      - ((1 + pinY[i] * gInv) * bpY[nID] - gInv * cpY[nID]) * apY[i] / bpY[nID] / bpY[nID];

    // Applying Net Weight
    float weight_this_net = netWeight[nID];
    pinGradX[i] *= weight_this_net;
    pinGradY[i] *= weight_this_net;
  }
}

__global__ void addPinGrad(
  const int    numPin,
  const int*   pin2Cell,
  const float* pinGradX,
  const float* pinGradY,
        float* cellGradX,
        float* cellGradY)
{
  const unsigned int pin_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(pin_id < numPin)
  {
    int cell_id = pin2Cell[pin_id];
    if(cell_id < 0)
      return;

    float pin_grad_x = pinGradX[pin_id];
    float pin_grad_y = pinGradY[pin_id];
    atomicAdd(&cellGradX[cell_id], pin_grad_x);
    atomicAdd(&cellGradY[cell_id], pin_grad_y);
  }
}

__global__ void updatePinCoordinateKernel(
  const int    numPin, 
  const int*   pin2Cell,
  const float* newCellX,
  const float* newCellY, 
  const float* pinOffsetX,
  const float* pinOffsetY,
        float* pinX,
        float* pinY)
{
  const unsigned int pin_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(pin_id < numPin)
  {
    int cell_id = pin2Cell[pin_id];

    // Some pins are from bterms (IO pins), 
    // so we should check if cell_id is valid.
    if(cell_id < 0)
      return;
    pinX[pin_id] = newCellX[cell_id] + pinOffsetX[pin_id];
    pinY[pin_id] = newCellY[cell_id] + pinOffsetY[pin_id];
  }
}

__global__ void computeNetBBox(
  const int    numNet,
  const int*   netStart,
  const float* pinX,
  const float* pinY,
        float* netBBoxWidth,
        float* netBBoxHeight)
{
  // i := netID
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < numNet)
  {
    float maxPinX = 0;
    float maxPinY = 0;
    float minPinX = k_float_max;
    float minPinY = k_float_max;

    // j = pin_id
    for(int j = netStart[i]; j < netStart[i+1]; j++)
    {
      maxPinX = max(pinX[j], maxPinX);
      minPinX = min(pinX[j], minPinX);
      maxPinY = max(pinY[j], maxPinY);
      minPinY = min(pinY[j], minPinY);
    }

    netBBoxWidth[i]  = maxPinX - minPinX;
    netBBoxHeight[i] = maxPinY - minPinY;
  }
}

__global__ void computeWAForEachNet(
  const int numNet,
  const float* bpX,
  const float* bpY,
  const float* bmX,
  const float* bmY,
  const float* cpX,
  const float* cpY,
  const float* cmX,
  const float* cmY,
        float* waForEachNetX,
        float* waForEachNetY)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < numNet)
  {
    waForEachNetX[i] = cpX[i] / bpX[i] - cmX[i] / bmX[i];
    waForEachNetY[i] = cpY[i] / bpY[i] - cmY[i] / bmY[i];
  }
}

WireLengthGradient::WireLengthGradient()
  : numCell_     (0),
    numNet_      (0),
    numPin_      (0),
    gammaInv_    (0.0),
    wl_grad_time_(0.0)
{}

WireLengthGradient::WireLengthGradient(std::shared_ptr<SkyPlaceDB> db)
  : WireLengthGradient()
{
  db_ = db;
  numCell_ = db_->numMovable();
  numNet_  = db_->numNet();
  numPin_  = db_->numPin();
  initForCUDAKernel();
}

void
WireLengthGradient::computeGrad(float* wlGradX, float* wlGradY)
{
  int numThread    = 64;
  int numBlockPin  = (numPin_      - 1 + numThread) / numThread; 
  int numBlockNet  = (numNet_  * 2 - 1 + numThread) / numThread; 

  auto wl_grad_start = getChronoNow();

  cudaDeviceSynchronize();

  d_bpX_.fillZero();
  d_bmX_.fillZero();
  d_bpY_.fillZero();
  d_bmY_.fillZero();

  d_cpX_.fillZero();
  d_cmX_.fillZero();
  d_cpY_.fillZero();
  d_cmY_.fillZero();

  cudaMemset(wlGradX, 0, sizeof(float) * numCell_);
  cudaMemset(wlGradY, 0, sizeof(float) * numCell_);

  // Step #1 : Compute Min Max
  computeMinMax<<<numBlockNet, numThread>>>(
    numNet_, 
    d_ptr_netStart_, 
    d_ptr_pinX_,
    d_ptr_pinY_,
    d_ptr_maxPinX_, 
    d_ptr_minPinX_, 
    d_ptr_maxPinY_, 
    d_ptr_minPinY_);

  // Step #2 : Compute a+ & a-, b+ & b-, c+ & c-
  computeABCKernel<<<numBlockPin, numThread>>>(
    numPin_, 
    gammaInv_, 
    d_ptr_pin2Net_, 
    d_ptr_pinX_, 
    d_ptr_pinY_, 
    d_ptr_maxPinX_, 
    d_ptr_maxPinY_, 
    d_ptr_minPinX_, 
    d_ptr_minPinY_,
    /* a+, a- */
    d_ptr_apX_, 
    d_ptr_apY_, 
    d_ptr_amX_, 
    d_ptr_amY_,
    /* b+, b- */
    d_ptr_bpX_, 
    d_ptr_bpY_, 
    d_ptr_bmX_, 
    d_ptr_bmY_,
    /* c+, c- */
    d_ptr_cpX_, 
    d_ptr_cpY_, 
    d_ptr_cmX_, 
    d_ptr_cmY_);

  // Step #3 : Compute Pin Gradient
  computePinGrad<<<numBlockPin, numThread>>>(
    numPin_,
    gammaInv_, 
    d_ptr_pin2Net_,
    d_ptr_pinX_, 
    d_ptr_pinY_,
    /* a+, a- */
    d_ptr_apX_, 
    d_ptr_apY_, 
    d_ptr_amX_, 
    d_ptr_amY_,
    /* b+, b- */
    d_ptr_bpX_, 
    d_ptr_bpY_, 
    d_ptr_bmX_, 
    d_ptr_bmY_,
    /* c+, c- */
    d_ptr_cpX_, 
    d_ptr_cpY_, 
    d_ptr_cmX_, 
    d_ptr_cmY_,
    d_ptr_netWeight_,
    d_ptr_pinGradX_, 
    d_ptr_pinGradY_);

  // Step #4 : Add each pin gradients of cells
  addPinGrad<<<numBlockPin, numThread>>>(
    numPin_, 
    d_ptr_pin2Cell_,
    d_ptr_pinGradX_, 
    d_ptr_pinGradY_, 
    wlGradX, 
    wlGradY);

  cudaDeviceSynchronize();
  const double time_this_iter = evalTime(wl_grad_start);
  wl_grad_time_ += time_this_iter;
}

void 
WireLengthGradient::updatePinCoordinates(const float* d_cellCx, const float* d_cellCy)
{
  int numThread = 64;
  int numBlockPin = (numPin_ - 1 + numThread) / numThread;

  updatePinCoordinateKernel<<<numBlockPin, numThread>>>(
    numPin_, 
    d_ptr_pin2Cell_,
    d_cellCx,
    d_cellCy,
    d_ptr_pinOffsetX_,
    d_ptr_pinOffsetY_,
    d_ptr_pinX_,
    d_ptr_pinY_);
}

float
WireLengthGradient::computeHPWL()
{
  int numThread = 64;
  int numBlockNet = (numNet_ - 1 + numThread) / numThread;

  computeNetBBox<<<numBlockNet, numThread>>>(
    numNet_, 
    d_ptr_netStart_,
    d_ptr_pinX_,
    d_ptr_pinY_,
    d_ptr_netBBoxWidth_, 
    d_ptr_netBBoxHeight_);

  float hpwl_x = computeVectorSum(d_netBBoxWidth_);
  float hpwl_y = computeVectorSum(d_netBBoxHeight_);

  return (hpwl_x + hpwl_y);
}

void
WireLengthGradient::initForCUDAKernel()
{
  // We don't want to store these data permanently
  std::vector<int> h_pin2Net(numPin_);
  std::vector<int> h_netStart(numNet_ + 1);
  std::vector<int> h_pin2Cell(numPin_, -1);

  std::vector<float> h_pinX(numPin_);
  std::vector<float> h_pinY(numPin_);

  std::vector<float> h_pinOffsetX(numPin_);
  std::vector<float> h_pinOffsetY(numPin_);

  std::vector<float> h_netWeight(numNet_);

  d_maxPinX_.resize(numNet_);
  d_minPinX_.resize(numNet_);

  d_maxPinY_.resize(numNet_);
  d_minPinY_.resize(numNet_);

  d_apX_.resize(numPin_); 
  d_bpX_.resize(numNet_);
  d_cpX_.resize(numNet_);

  d_apY_.resize(numPin_);
  d_bpY_.resize(numNet_);
  d_cpY_.resize(numNet_);

  d_amX_.resize(numPin_);
  d_bmX_.resize(numNet_);
  d_cmX_.resize(numNet_);

  d_amY_.resize(numPin_);
  d_bmY_.resize(numNet_);
  d_cmY_.resize(numNet_);

  d_pinGradX_.resize(numPin_);
  d_pinGradY_.resize(numPin_);

  d_netBBoxWidth_ .resize(numNet_);
  d_netBBoxHeight_.resize(numNet_);

  d_waForEachNetX_.resize(numNet_);
  d_waForEachNetY_.resize(numNet_);

  d_pin2Net_.resize(numPin_);
  d_pin2Cell_.resize(numPin_);
  d_netStart_.resize(numNet_ + 1);

  d_pinX_.resize(numPin_);
  d_pinY_.resize(numPin_);

  d_pinOffsetX_.resize(numPin_);
  d_pinOffsetY_.resize(numPin_);

  d_netWeight_.resize(numNet_);

  d_ptr_pinX_ = d_pinX_.data();
  d_ptr_pinY_ = d_pinY_.data();
  
  d_ptr_pin2Net_        = d_pin2Net_.data();
  d_ptr_pin2Cell_       = d_pin2Cell_.data();
  d_ptr_netStart_       = d_netStart_.data();

  d_ptr_pinOffsetX_     = d_pinOffsetX_.data();
  d_ptr_pinOffsetY_     = d_pinOffsetY_.data();
  d_ptr_netBBoxWidth_   = d_netBBoxWidth_.data();
  d_ptr_netBBoxHeight_  = d_netBBoxHeight_.data();
  d_ptr_waForEachNetX_  = d_waForEachNetX_.data();
  d_ptr_waForEachNetY_  = d_waForEachNetY_.data();

  d_ptr_maxPinX_        = d_maxPinX_.data();
  d_ptr_minPinX_        = d_minPinX_.data();

  d_ptr_maxPinY_        = d_maxPinY_.data();
  d_ptr_minPinY_        = d_minPinY_.data();

  d_ptr_apX_            = d_apX_.data();
  d_ptr_bpX_            = d_bpX_.data();
  d_ptr_cpX_            = d_cpX_.data();

  d_ptr_apY_            = d_apY_.data();
  d_ptr_bpY_            = d_bpY_.data();
  d_ptr_cpY_            = d_cpY_.data();

  d_ptr_amX_            = d_amX_.data();
  d_ptr_bmX_            = d_bmX_.data();
  d_ptr_cmX_            = d_cmX_.data();

  d_ptr_amY_            = d_amY_.data();
  d_ptr_bmY_            = d_bmY_.data();
  d_ptr_cmY_            = d_cmY_.data();

  d_ptr_pinGradX_       = d_pinGradX_.data();
  d_ptr_pinGradY_       = d_pinGradY_.data();

  d_ptr_netWeight_      = d_netWeight_.data();

  // pinListCell is a list of pin,
  // whose order is determined by cells
  // [Example]
  //   pinListCell[0] : pin_id of 1st pin of cell 0
  //   pinListCell[1] : pin_id of 2nd pin of cell 0
  //   pinListCell[2] : pin_id of 3rd pin of cell 0

  //   pinListCell[3] : pin_id of 1st pin of cell 1
  //   pinListCell[4] : pin_id of 2nd pin of cell 1
  //   pinListCell[5] : pin_id of 3rd pin of cell 1

  //   ...

  // cellNumPin contains the number of pins 
  // in each cells as a list
  // Since the order of PinID is aligned with NetID 
  // We need extra data structure so that we can access 
  // to the correct PinID with Cell ID

  // Build Netlist Information Vectors
  for(auto& net : db_->nets())
  {
    int netID = net->id();
    h_netStart[netID]  = net->pins()[0]->id();
    h_netWeight[netID] = net->weight();
    for(auto& pin : net->pins())
    {
      int pin_id = pin->id();
      h_pin2Net[pin_id] = netID;
    }
  }

  h_netStart[numNet_] = numPin_;
  // we should not leave netStart_[numNet_ + 1]
  // as a garbage value

  for(const auto& cell : db_->movableCells())
  {
    int cell_id = cell->id();
    for(auto& pin : cell->pins())
    {
      int pin_id = pin->id();
      h_pin2Cell[pin_id] = cell_id;
    }
  }

  for(auto& pin : db_->pins())
  {
    int pin_id = pin->id();
    h_pinX[pin_id] = pin->cx();
    h_pinY[pin_id] = pin->cy();
    h_pinOffsetX[pin_id] = pin->offsetX();
    h_pinOffsetY[pin_id] = pin->offsetY();
  }

  d_pinX_        = h_pinX;
  d_pinY_        = h_pinY;
  d_pin2Net_     = h_pin2Net;
  d_pin2Cell_    = h_pin2Cell;
  d_netStart_    = h_netStart;
  d_pinOffsetX_  = h_pinOffsetX;
  d_pinOffsetY_  = h_pinOffsetY;
  d_netWeight_   = h_netWeight;
}

}; // namespace skyplace
