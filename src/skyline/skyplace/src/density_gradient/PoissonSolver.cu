#include "PoissonSolver.h"
#include "CudaDct.h"

#define FFT_PI 3.141592653589793238462643383279502884197169

namespace skyplace 
{

__global__ void precomputeExpk(cufftComplex *expkM, 
                               cufftComplex *expkN,
                               const int M, 
                               const int N)
{
  const int tID = blockDim.x * blockIdx.x + threadIdx.x;
  
  if(tID <= M / 2)
  {
    int hID = tID;
    cufftComplex W_h_4M = make_float2(__cosf((float)FFT_PI * hID / (2 * M)), 
                                     -__sinf((float)FFT_PI * hID / (M * 2)));
    expkM[hID] = W_h_4M;
  }
  if(tID <= N / 2)
  {
    int wid = tID;
    cufftComplex W_w_4N = make_float2(__cosf((float)FFT_PI * wid / (2 * N)), 
                                     -__sinf((float)FFT_PI * wid / (N * 2)));
    expkN[wid] = W_w_4N;
  }
}

__global__ void precomputeExpkForInverse(cufftComplex *expkM, 
                                         cufftComplex *expkN, 
                                         cufftComplex *expkMN_1, 
                                         cufftComplex *expkMN_2, 
                                         const int M, 
                                         const int N)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if(tid < M)
  {
    int hid = tid;
    cufftComplex W_h_4M 
      = make_float2(__cosf((float)FFT_PI * hid / (2 * M)), 
                  - __sinf((float)FFT_PI * hid / (M * 2)));
    expkM[hid] = W_h_4M;
    // expkMN_1
    cufftComplex W_h_4M_offset 
      = make_float2(__cosf((float)FFT_PI * (hid + M) / (2 * M)), 
                  - __sinf((float)FFT_PI * (hid + M) / (M * 2)));
    expkMN_1[hid] = W_h_4M;
    expkMN_1[hid + M] = W_h_4M_offset;
  
    //expkMN_2
    W_h_4M = make_float2(-__sinf((float)FFT_PI * (hid - (N - 1)) / (M * 2)), 
                         -__cosf((float)FFT_PI * (hid - (N - 1)) / (2 * M)));

    W_h_4M_offset = make_float2(-__sinf((float)FFT_PI * (hid - (N - 1) + M) / (M * 2)), 
                                -__cosf((float)FFT_PI * (hid - (N - 1) + M) / (2 * M)));

    expkMN_2[hid] = W_h_4M;
    expkMN_2[hid + M] = W_h_4M_offset;
  }
  if(tid <= N / 2)
  {
    int wid = tid;
    cufftComplex W_w_4N = make_float2(__cosf((float)FFT_PI * wid / (2 * N)), 
                                     -__sinf((float)FFT_PI * wid / (N * 2)));
    expkN[wid] = W_w_4N;
  }
}

__global__ void divideByWSquare(const int numBinX,
                                const int numBinY,
                                cufftReal* input)
{
  //const int binID = blockDim.x * blockIdx.x + threadIdx.x;

  const int wID = blockDim.x * blockIdx.x + threadIdx.x;
  const int hID = blockDim.y * blockIdx.y + threadIdx.y;

  if(wID < numBinX && hID < numBinY)
  {
    int binID = wID + hID * numBinX;

    if(hID == 0 && wID == 0)
      input[binID] = 0.0;
    else
    {
      float denom1 = (2.0 * float(FFT_PI) * wID) / numBinX; 
      float denom2 = (2.0 * float(FFT_PI) * hID) / numBinY;
      input[binID] /= (denom1 * denom1 + denom2 * denom2);
    }
  }
}

__global__ void multiplyW(
  const int numBinX,
  const int numBinY,
  const cufftReal* auv, 
        cufftReal* inputForX,
        cufftReal* inputForY)
{
  //const int binID = blockDim.x * blockIdx.x + threadIdx.x;
  const int wID = blockDim.x * blockIdx.x + threadIdx.x;
  const int hID = blockDim.y * blockIdx.y + threadIdx.y;

  if(wID < numBinX && hID < numBinY)
  {
    int binID = wID + hID * numBinX;

    float w_u = (2.0 * float(FFT_PI) * wID) / numBinX; 
    float w_v = (2.0 * float(FFT_PI) * hID) / numBinY; 

    inputForX[binID] = w_u * auv[binID]; // why w_v?
    inputForY[binID] = w_v * auv[binID]; // why w_u?
  }
}

PoissonSolver::PoissonSolver() : numBinX_(0), numBinY_(0) {}

PoissonSolver::PoissonSolver(int numBinX, int numBinY) : PoissonSolver()
{
  numBinX_ = numBinX;
  numBinY_ = numBinY;
  init();
}

void
PoissonSolver::init()
{
  setupForCUDAKernel();
}

void
PoissonSolver::solvePoissonPotential(const float* binDensity, float* potential)
{
  int numThread = 16;

  dim3 gridSize( (numBinX_ + numThread - 1) / numThread, 
                 (numBinY_ + numThread - 1) / numThread, 1);

  dim3 blockSize(numThread, numThread, 1);

  // Step #1. Compute Coefficient (a_uv)
  dct_2d_fft(numBinY_, 
             numBinX_,
             plan_,
             d_expkM_,
             d_expkN_,
             binDensity, 
             d_workSpaceReal1_, 
             d_workSpaceComplex_, 
             d_auv_);

  // Step #2. Divide by (w_u^2 + w_v^2)
  divideByWSquare<<<gridSize, blockSize>>>(numBinX_, numBinY_, d_auv_);

  // Step #3. Compute Potential
  idct_2d_fft(
    numBinY_, 
    numBinX_,
    planInverse_,
    d_expkMForInverse_,
    d_expkNForInverse_,
    d_expkMN1_,
    d_expkMN2_,
    d_auv_, 
    d_workSpaceComplex_, 
    d_workSpaceReal1_, 
    potential);
}

void
PoissonSolver::solvePoissonForce(
  float* electroForceX, 
  float* electroForceY)
{
  int numThread = 16;

  dim3 gridSize( (numBinX_ + numThread - 1) / numThread, 
                 (numBinY_ + numThread - 1) / numThread, 1);

  dim3 blockSize(numThread, numThread, 1);

  // Step #4. Multiply w_u , w_v
  multiplyW<<<gridSize, blockSize>>>(numBinX_, 
                                     numBinY_,
                                     d_auv_, 
                                     d_inputForX_, 
                                     d_inputForY_);

  // Step #5. Compute ElectroForceX
  idxst_idct(numBinX_, 
             numBinY_,
             planInverse_,
             d_expkMForInverse_,
             d_expkNForInverse_,
             d_expkMN1_,
             d_expkMN2_,
             d_inputForX_, 
             d_workSpaceReal1_, 
             d_workSpaceComplex_, 
             d_workSpaceReal2_, 
             d_workSpaceReal3_, 
             electroForceX);

  // Step #6. Compute ElectroForceY
  idct_idxst(numBinY_, 
             numBinX_,
             planInverse_,
             d_expkMForInverse_,
             d_expkNForInverse_,
             d_expkMN1_,
             d_expkMN2_,
             d_inputForY_, 
             d_workSpaceReal1_, 
              d_workSpaceComplex_, 
             d_workSpaceReal2_, 
             d_workSpaceReal3_, 
             electroForceY);
}

void
PoissonSolver::solvePoisson(const float* binDensity,
                                  float* potential,
                                  float* electroForceX,
                                  float* electroForceY)
{
  int numThread = 16;

  dim3 gridSize( (numBinX_ + numThread - 1) / numThread, 
                 (numBinY_ + numThread - 1) / numThread, 1);

  dim3 blockSize(numThread, numThread, 1);

  // Step #1. Compute Coefficient (a_uv)
  dct_2d_fft(numBinY_, 
             numBinX_,
             plan_,
             d_expkM_,
             d_expkN_,
             binDensity, 
             d_workSpaceReal1_, 
             d_workSpaceComplex_, 
             d_auv_);

  // Step #2. Divide by (w_u^2 + w_v^2)
  divideByWSquare<<<gridSize, blockSize>>>(numBinX_, numBinY_, d_auv_);

  // Step #3. Compute Potential
  idct_2d_fft(numBinY_, 
              numBinX_,
              planInverse_,
              d_expkMForInverse_,
              d_expkNForInverse_,
              d_expkMN1_,
              d_expkMN2_,
              d_auv_, 
              d_workSpaceComplex_, 
              d_workSpaceReal1_, 
              potential);

  // Step #4. Multiply w_u , w_v
  multiplyW<<<gridSize, blockSize>>>(numBinX_, 
                                     numBinY_,
                                     d_auv_, 
                                     d_inputForX_, 
                                     d_inputForY_);

  // Step #5. Compute ElectroForceX
  idxst_idct(numBinY_, 
             numBinX_,
             planInverse_,
             d_expkMForInverse_,
             d_expkNForInverse_,
             d_expkMN1_,
             d_expkMN2_,
             d_inputForX_, 
             d_workSpaceReal1_, 
             d_workSpaceComplex_, 
             d_workSpaceReal2_, 
             d_workSpaceReal3_, 
             electroForceX);

  // Step #6. Compute ElectroForceY
  idct_idxst(numBinY_, 
             numBinX_,
             planInverse_,
             d_expkMForInverse_,
             d_expkNForInverse_,
             d_expkMN1_,
             d_expkMN2_,
             d_inputForY_, 
             d_workSpaceReal1_, 
             d_workSpaceComplex_, 
             d_workSpaceReal2_, 
             d_workSpaceReal3_, 
             electroForceY);

  cudaDeviceSynchronize();
}

void
PoissonSolver::setupForCUDAKernel()
{
  const int num_bin_2 = numBinX_ * numBinY_;

  cudaMalloc((void**)&d_binDensity_, num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_auv_,        num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_potential_,  num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_efX_,        num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_efY_,        num_bin_2 * sizeof(cufftReal));

  cudaMalloc((void**)&d_workSpaceReal1_, num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_workSpaceReal2_, num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_workSpaceReal3_, num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void**)&d_workSpaceComplex_, (numBinX_ / 2 + 1) * numBinY_ * sizeof(cufftComplex));

  // expk
  // For DCT2D
  cudaMalloc((void **)&d_expkM_, (numBinY_ / 2 + 1) * sizeof(cufftComplex));
  cudaMalloc((void **)&d_expkN_, (numBinX_ / 2 + 1) * sizeof(cufftComplex));

  // For IDCT2D & IDXST_IDCT & IDCT_IDXST
  cudaMalloc((void **)&d_expkMForInverse_, (numBinY_) * sizeof(cufftComplex));
  cudaMalloc((void **)&d_expkNForInverse_, (numBinX_ / 2 + 1) * sizeof(cufftComplex));
  cudaMalloc((void **)&d_expkMN1_, (numBinX_ + numBinY_) * sizeof(cufftComplex));
  cudaMalloc((void **)&d_expkMN2_, (numBinX_ + numBinY_) * sizeof(cufftComplex));

  // For Input For IDXST_IDCT & IDCT_IDXST
  cudaMalloc((void **)&d_inputForX_, num_bin_2 * sizeof(cufftReal));
  cudaMalloc((void **)&d_inputForY_, num_bin_2 * sizeof(cufftReal));

  int numThread = 1024;
  int numBin = std::max(numBinX_, numBinY_);
  int numBlock  = (numBin - 1 + numThread) / numThread;

  precomputeExpk<<<numBlock, numThread>>>(d_expkM_, 
                                          d_expkN_, 
                                          numBinY_, 
                                          numBinX_);

  precomputeExpkForInverse<<<numBlock, numThread>>>(d_expkMForInverse_, 
                                                    d_expkNForInverse_,
                                                    d_expkMN1_, 
                                                    d_expkMN2_,
                                                    numBinY_, 
                                                    numBinX_);

  cufftPlan2d(&plan_,        numBinY_, numBinX_, CUFFT_R2C);
  cufftPlan2d(&planInverse_, numBinY_, numBinX_, CUFFT_C2R);
}

void
PoissonSolver::freeDeviceMemory()
{
  cudaFree(d_binDensity_);
  cudaFree(d_auv_);
  cudaFree(d_potential_);

  cudaFree(d_efX_);
  cudaFree(d_efY_);

  cudaFree(d_workSpaceReal1_);
  cudaFree(d_workSpaceReal2_);
  cudaFree(d_workSpaceReal3_);

  cudaFree(d_workSpaceComplex_);

  cudaFree(d_expkN_);
  cudaFree(d_expkM_);

  cudaFree(d_expkNForInverse_);
  cudaFree(d_expkMForInverse_);

  cudaFree(d_expkMN1_);
  cudaFree(d_expkMN2_);

  cudaFree(d_inputForX_);
  cudaFree(d_inputForY_);

  cufftDestroy(plan_);
  cufftDestroy(planInverse_);
}

} // namespace skyplace
