#ifndef INITIAL_PLACE_H
#define INITIAL_PLACE_H

#include <memory>

#include "cuda_linalg/EigenDef.h"
#include "gpu_louvain/gpulouvain.h"

#include "HyperParam.h" // InitialPlaceMethod

namespace skyplace 
{

using namespace cuda_linalg;
using namespace Louvain;

class GPCell;
class SkyPlaceDB;

struct InitialPlaceParams
{
  int num_cpu_threads;
  int sparse_threshold;
  int max_degree_cluster;
  int max_degree_laplacian;
  double alpha_x;
  double alpha_y;
  double filler_dev_coeff_x;
  double filler_dev_coeff_y;

  double random_init_dev_coeff_x;
  double random_init_dev_coeff_y;
};

class InitialPlacer
{
  public:
    
    InitialPlacer();                               // Default Constructor
    InitialPlacer(std::shared_ptr<SkyPlaceDB> db); // Constructor really used

    // APIs
    void doInitialPlace(InitialPlaceMethod method);

    // Getters
    double getRuntime() const { return runtime_; }

  private:

    std::shared_ptr<SkyPlaceDB> db_;
    InitialPlaceParams ip_param_;

    EigenVector  xm_;     // Movable x
    EigenVector  ym_;     // Movalbe y

    EigenVector  xf_;     // Fixed x
    EigenVector  yf_;     // Fixed y

    EigenVector  Lmf_xf_; // Lmf * xf
    EigenVector  Lmf_yf_; // Lmf * yf

    EigenSMatrix L_;      // Full Laplacian
    EigenSMatrix Lmm_;    // Laplacian between movable cells
    EigenSMatrix Lff_;    // Laplacian between fixed   cells

    void doRandomInit();
    void doSdpInit();

    double runtime_;
    double total_cluster_area_;

    // Clustering-based Initialization
    int num_cluster_;
    std::vector<double> v_vector_;

    std::unordered_map<int, int> after2before_;
    std::unordered_map<int, int> before2after_;

    int garbageClusterID_; // To ignore sparse cluster

    // Step #1
    HostStructure host_structure_;
    void buildHostStructure();
    void doClustering(HostStructure& host_structure);
    void refineCluster();
    void ignoreSparse(int minSize);
    void computeSizeMap();
    void computeFixedInfo();
    void simpleBiPartitioning(std::vector<GPCell*>& smallGraph, int& numCluster);

    // Step #2
    void createClusterLaplacian(EigenSMatrix& L);
    void extractPartialLaplacian(
      const EigenVector&  xf,
      const EigenVector&  yf,
      const EigenSMatrix& L, 
            EigenSMatrix& Lff,
            EigenSMatrix& Lmm,
            EigenVector&  Lmf_xf,
            EigenVector&  Lmf_yf);

    // Step #3
    EigenVector solveSDP(
      const EigenSMatrix& Lmm,
      const EigenVector& Lmf_xf, // or Lmf_yf
      const std::vector<double>& v_vector, 
      const double K);

    double DbToSdpX(double llx) const;
    double DbToSdpY(double lly) const;

    double SdpToDbX(double ccx) const;
    double SdpToDbY(double ccy) const;

    double getMirrorX(double locX, double dieCx, double dieLx, double dieUx) const;
    double getMirrorY(double locY, double dieCy, double dieLy, double dieUy) const;
};

}; // namespace skyplace 

#endif
