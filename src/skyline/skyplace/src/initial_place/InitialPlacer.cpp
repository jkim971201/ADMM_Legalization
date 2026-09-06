#include <cstdio>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <memory>
#include <chrono>

#include "SkyPlaceDB.h"
#include "InitialPlacer.h"
#include "object/GPObject.h"
#include "Util.h"

#include "SDPSolverCPU.h"
//#include "SDPSolverGPU.h"

#include "FMPartitioner.h"

namespace skyplace 
{

inline bool isLargeMacro(GPCell* cell, double dieDx, double dieDy)
{
  return (cell->dx() > dieDx * 0.2) || (cell->dy() > dieDy * 0.2);
}

InitialPlacer::InitialPlacer()
  : db_         (nullptr),
    num_cluster_(0),
    runtime_    (0)
{}

InitialPlacer::InitialPlacer(std::shared_ptr<SkyPlaceDB> db)
  : InitialPlacer()
{
  db_ = db;
  
  // Initialize Hyperparameters
  ip_param_.num_cpu_threads = 4; // 8 causes crash?

  ip_param_.sparse_threshold     = 3;
  ip_param_.max_degree_cluster   = 24;
  ip_param_.max_degree_laplacian = 25;
  ip_param_.alpha_x = 10000000;
  ip_param_.alpha_y = 10000000;
  ip_param_.filler_dev_coeff_x = 0.3;
  ip_param_.filler_dev_coeff_y = 0.3;

  ip_param_.random_init_dev_coeff_x = 0.05;
  ip_param_.random_init_dev_coeff_y = 0.05;
}

void
InitialPlacer::doInitialPlace(InitialPlaceMethod method)
{
  auto initial_place_start = std::chrono::high_resolution_clock::now();

  const std::string ip_method
    = method == InitialPlaceMethod::RANDOM_START ? "random_start"
                                                 : "sdp_relaxation";

  printf("\n");
  printf("Run initial_place... (%s)\n", ip_method.c_str());

  if(method == InitialPlaceMethod::RANDOM_START)
    doRandomInit();
  else 
    doSdpInit();

  db_->updateHpwl();

  runtime_ = evalTime(initial_place_start);
  printf("Initial HPWL: %.1f\n", db_->getHPWL() / 1e6);
  printf("--> InitialPlacement takes %5.2f s\n", runtime_);
}

void
InitialPlacer::doRandomInit()
{
  // Place all movable cells in the center with Gaussian Noise.
  // See DREAMPlace TCAD`21 for details
  // The authors argue that random center initiailization
  // does not degrade the placement quality
  const double die_cx = db_->die()->cx();
  const double die_cy = db_->die()->cy();

  const double die_lx = db_->die()->lx();
  const double die_ly = db_->die()->ly();

  const double die_ux = db_->die()->ux();
  const double die_uy = db_->die()->uy();

  const double die_dx = db_->die()->dx();
  const double die_dy = db_->die()->dy();

  double mean_x = die_cx;  
  double mean_y = die_cy;  

  double deviation_x = die_dx * ip_param_.random_init_dev_coeff_x; // 0.05
  double deviation_y = die_dy * ip_param_.random_init_dev_coeff_y; // 0.05

  std::default_random_engine gen;
  std::normal_distribution<double> filler_noise_x(mean_x, deviation_x);
  std::normal_distribution<double> filler_noise_y(mean_y, deviation_y);

  auto& movable_cells = db_->movableCells();
  const int num_movable = movable_cells.size();

  #pragma omp parallel for num_threads(ip_param_.num_cpu_threads)
  for(int i = 0; i < num_movable; i++)
  {
    GPCell* cell = movable_cells[i];
    double loc_x = filler_noise_x(gen);
    double loc_y = filler_noise_y(gen);
    if(cell->isFiller() == true)
    {
      loc_x = getMirrorX(loc_x, die_cx, die_lx, die_ux);
      loc_y = getMirrorY(loc_y, die_cy, die_ly, die_uy);
    }
    cell->setCenterLocation(loc_x, loc_y);
    db_->moveCellInsideLayout(cell);
  }
}

void 
InitialPlacer::simpleBiPartitioning(std::vector<GPCell*>& smallGraph, int& numCluster)
{
  FMPartitioner fm(smallGraph, 0.3, numCluster);
}

void 
InitialPlacer::refineCluster()
{
  std::unordered_map<int, std::vector<GPCell*>> clusterIDtoCells;

  double dieDx = db_->die()->dx();
  double dieDy = db_->die()->dy();

  const auto& movable_cells = db_->movableCells();
  for(auto& cell : movable_cells)
  {
    if(!cell->isFiller() && cell->pins().size() != 0)
    {
      int original_cluster_id = cell->clusterID();
      clusterIDtoCells[original_cluster_id].push_back(cell);
    }
  }

  for(auto& [cluster_id, cells] : clusterIDtoCells)
  {
    int numMacro = 0;
    int numCell  = cells.size();

    for(auto& cell : cells)
    {
      if(cell->isMacro() == true)
        numMacro++;
      if(isLargeMacro(cell, dieDx, dieDy) && numCell != 1)
      {
        cell->setClusterID(num_cluster_);
        num_cluster_++;
      }
    }

    if(numMacro > 15 && numMacro < 100)
    {
      // printf("Too many macros are detected in a cluster...\n");
      std::vector<GPCell*> smallGraph; // only has Macros
      
      for(auto& cell : cells)
      {
        if(cell->isMacro())
          smallGraph.push_back(cell);
      }
      simpleBiPartitioning(smallGraph, num_cluster_);
    }
  }

  ignoreSparse(ip_param_.sparse_threshold);
}

void
InitialPlacer::buildHostStructure() 
{
  const int V = db_->numMovable() - db_->numFiller();

  host_structure_.originalV = V;
  host_structure_.V = V;

  cudaHostAlloc((void**)&host_structure_.vertexCommunity, V * sizeof(int), cudaHostAllocDefault);
  cudaHostAlloc((void**)&host_structure_.communityWeight, V * sizeof(float), cudaHostAllocDefault);
  cudaHostAlloc((void**)&host_structure_.edgesIndex, (V + 1) * sizeof(int), cudaHostAllocDefault);
  cudaHostAlloc((void**)&host_structure_.originalToCommunity, V * sizeof(int), cudaHostAllocDefault);

  std::vector<std::unordered_map<int, float>> neighbours(V);
  // Here is assumption that graph is undirected.

  float* host_community_weight = host_structure_.communityWeight;

  const auto& gp_nets = db_->nets();
  const int num_nets = gp_nets.size();
  const int max_degree = ip_param_.max_degree_cluster;

  #pragma omp parallel for num_threads(ip_param_.num_cpu_threads)
  for(int net_id = 0; net_id < num_nets; net_id++)
  {
    GPNet* net = gp_nets[net_id];
    const int net_degree = net->deg();
    if(net_degree < 2 || net_degree > max_degree) // 25?
      continue;

    float clique_weight = 1.0 / (float(net->deg()) - 1.0);
    const auto& net_pins = net->pins();
    for(int p1 = 0; p1 < net_degree - 1; p1++)
    {
      GPPin*   pin1 = net_pins[p1];
      GPCell* cell1 = pin1->cell(); 

      if(pin1->isIO() || cell1->isFixed())
        continue;

      int v1 = cell1->id();
      for(int p2 = p1 + 1; p2 < net_degree; p2++)
      {
        GPPin*   pin2 = net_pins[p2];
        GPCell* cell2 = pin2->cell(); 
        
        int v2 = cell2->id();
        if(v1 == v2 || pin2->isIO() || cell2->isFixed())
          continue;

        #pragma omp atomic
        host_community_weight[v1] += clique_weight;

        #pragma omp atomic
        host_community_weight[v2] += clique_weight;

        #pragma omp atomic
        host_structure_.M += clique_weight;

        #pragma omp atomic
        neighbours[v1][v2] += clique_weight;

        #pragma omp atomic
        neighbours[v2][v1] += clique_weight;
      }
    }
  }

  int E = 0;
  for(const auto& m : neighbours)
    E += m.size();

  host_structure_.E = E;
  cudaHostAlloc((void**)&host_structure_.edges,   E * sizeof(  int), cudaHostAllocDefault);
  cudaHostAlloc((void**)&host_structure_.weights, E * sizeof(float), cudaHostAllocDefault);

  int edge_index = 0;
  for(int vertex_id = 0; vertex_id < V; vertex_id++) 
  {
    host_structure_.edgesIndex[vertex_id] = edge_index;
    const auto& pairs_this_vertex = neighbours[vertex_id];
    for(const auto& [opposite_vertex_id, weight] : pairs_this_vertex)
    {
      host_structure_.edges[edge_index] = opposite_vertex_id;
      host_structure_.weights[edge_index] = weight;
      edge_index++;
    }
  }

  host_structure_.edgesIndex[V] = E;
  printf("Original Graph (V, E) = (%d, %d)\n", V, E);
}

void
InitialPlacer::doClustering(HostStructure& host_structure)
{
  // GPU Louvain Clustering
  gpuLouvain(host_structure, /* min_gain */ 0.01, ip_param_.sparse_threshold /* min_cluster_size */);
  num_cluster_ = host_structure.num_vertex_final;

  for(auto& cell : db_->movableCells() )
  {
    if(cell->isFiller() == true)
      continue;
    assert(cell->pins().size() != 0);
    int cluster_id = host_structure_.originalToCommunity[cell->id()];
    cell->setClusterID(cluster_id);
  }
}

void
InitialPlacer::computeFixedInfo()
{
  const int num_fixed_inst = db_->numFixed();
  const int num_io         = db_->numIO();
  const int num_fixed      = num_fixed_inst + num_io;

  xf_.resize(num_fixed);
  yf_.resize(num_fixed);

  int fixed_id = 0;
  for(const auto& cell : db_->fixedCells())
  {
    xf_(fixed_id) = DbToSdpX(cell->cx());
    yf_(fixed_id) = DbToSdpY(cell->cy());
    fixed_id++;
  }

  for(const auto& io : db_->getIOPins())
  {
    xf_(fixed_id) = DbToSdpX(io->cx());
    yf_(fixed_id) = DbToSdpY(io->cy());
    fixed_id++;
  }
}

void
InitialPlacer::createClusterLaplacian(EigenSMatrix& L)
{
  std::unordered_map<GPPin*, int> pin_to_vertex_id;
  const int num_fixed_inst = db_->numFixed();

  int io_id = num_cluster_ + num_fixed_inst;
  for(auto& io_pin_ptr : db_->getIOPins())
    pin_to_vertex_id[io_pin_ptr] = io_id++;

  // vertex_id - verteix_id - edge_weight
  std::vector<EigenTriplet> triplet_vector;

  const int num_total_vertex = num_cluster_ + xf_.size();
  // M2M : Movable to Movable
  // M2F : Movable to Fixed
  constexpr double k_weight_for_M2M =  1.0;
  constexpr double k_weight_for_M2F = 50.0;

  const int max_degree = ip_param_.max_degree_laplacian;

  auto loop_start = std::chrono::high_resolution_clock::now();
  for(const auto& net : db_->nets())
  {
    const int net_degree = net->deg();
    // Since we are ignoring too large nets,
    // there can be some clusters that are not 
    // connected to any other vertex.
    if(net_degree < 2 || net_degree > max_degree)
      continue;

    const auto& net_pins = net->pins();
    for(int p1 = 0; p1 < net_degree - 1; p1++)
    {
      GPPin* pin1 = net_pins[p1];

      int v1 = -1;
      bool v1_fixed = false;
      if(pin1->isIO() == true)
      {
        auto find_vertex_id = pin_to_vertex_id.find(pin1);
        assert(find_vertex_id != pin_to_vertex_id.end());
        v1 = find_vertex_id->second;
        v1_fixed = true;
      }
      else
      {
        GPCell* cell1 = pin1->cell(); 
        if(cell1->isFixed() == true)
        {
          v1 = cell1->id() + num_cluster_;
          v1_fixed = true;
        }
        else
          v1 = cell1->clusterID();
      }
      
      for(int p2 = p1 + 1; p2 < net_degree; p2++)
      {
        GPPin* pin2 = net_pins[p2];

        int v2 = -1;
        bool v2_fixed = false;
        if(pin2->isIO() == true)
        {
          auto find_vertex_id = pin_to_vertex_id.find(pin2);
          assert(find_vertex_id != pin_to_vertex_id.end());
          v2 = find_vertex_id->second;
          v2_fixed = true;
        }
        else
        {
          GPCell* cell2 = pin2->cell(); 
          if(cell2->isFixed() == true)
          {
            v2 = cell2->id() + num_cluster_;
            v2_fixed = true;
          }
          else
            v2 = cell2->clusterID();
        }
        
        if(v1 == v2)
          continue;

        double weight = 1.0;
        if(v1_fixed == true || v2_fixed == true)
          weight = k_weight_for_M2F;
        else
          weight = k_weight_for_M2M;

        // L = D - A
        // For Diagonal Matrix
        triplet_vector.push_back(EigenTriplet(v1, v1, +weight));
        triplet_vector.push_back(EigenTriplet(v2, v2, +weight));

        // For Adjacency Matrix
        triplet_vector.push_back(EigenTriplet(v1, v2, -weight));
        triplet_vector.push_back(EigenTriplet(v2, v1, -weight));
      }
    }
  }

  L.resize(num_total_vertex, num_total_vertex);
  L.setFromTriplets(triplet_vector.begin(), triplet_vector.end());
}

void
InitialPlacer::extractPartialLaplacian(
  const EigenVector&  xf,
  const EigenVector&  yf,
  const EigenSMatrix& L, // Full Laplacian
        EigenSMatrix& Lff,
        EigenSMatrix& Lmm,
        EigenVector&  Lmf_xf,
        EigenVector&  Lmf_yf)
{
  const int num_fixed = xf.size();

  EigenSMatrix Lmf 
    = L.block(0,            // Start Index of Row (Y)
              num_cluster_, // Start Index of Col (X)
              num_cluster_, // Size of Y (numRow)
              num_fixed);   // Size of X (numCol)

  Lmm = L.block(0, 0, num_cluster_, num_cluster_);
  Lff = L.block(num_cluster_, num_cluster_, num_fixed, num_fixed);

  Lmf_xf.resize(num_cluster_);
  Lmf_yf.resize(num_cluster_);

  Lmf_xf = Lmf * xf;
  Lmf_yf = Lmf * yf;
}

EigenVector 
InitialPlacer::solveSDP(
  const EigenSMatrix&        Lmm, 
  const EigenVector&         b, 
  const std::vector<double>& v, 
  const double               K) 
{
  int N = Lmm.rows();

  auto t1 = getChronoNow();
  EigenSMatrix M_0;           // Objective
  EigenSMatrix M_1, M_2, M_3; // Constraint
  M_0.resize(N + 1, N + 1);
  M_1.resize(N + 1, N + 1);
  M_2.resize(N + 1, N + 1);
  M_3.resize(N + 1, N + 1);

  for(int i = 0; i < N; i++)
  {
    M_0.coeffRef(i + 1, i + 1) = Lmm.coeff(i, i);
    M_0.coeffRef(i + 1, 0) = b(i);
    M_0.coeffRef(0, i + 1) = b(i);

    M_1.coeffRef(i + 1, i + 1) = v[i];

    M_2.coeffRef(0, i + 1) = v[i]; 
    M_2.coeffRef(i + 1, 0) = v[i]; 
  }
  
  // M_3 has only one non-zero in (0, 0)
  M_3.coeffRef(0, 0) = 1;

  sdp_solver::SDPSolverCPU solver(N + 1);
  //sdp_solver::SDPSolverGPU solver(N + 1);
  solver.setObjective(M_0);
  solver.addEqualityConstraint(M_1, K);
  solver.addEqualityConstraint(M_2, 0);
  solver.addEqualityConstraint(M_3, 1);
  solver.solve();

  printf("ObjVal : %f\n", solver.getObjectiveValue());

  return solver.getResult();
}

void 
InitialPlacer::computeSizeMap()
{
  v_vector_.resize(num_cluster_, 0.0);
  for(auto& cell : db_->movableCells() )
  {
    if(!cell->isFiller())
    {
      int cluster_id = cell->clusterID();
      const double width  = static_cast<double>(cell->dx());
      const double height = static_cast<double>(cell->dy());
      const double area = width * height;
      v_vector_[cluster_id] += area;
      total_cluster_area_ += area;
    }
  }

//  double avg_cluster_area 
//    = total_cluster_area_ / num_cluster_;
//  for(int i = 0; i < num_cluster_; i++)
//  {
//    double area_this_cluster = v_vector_[i];
//    v_vector_[i] = area_this_cluster / avg_cluster_area;
//  }
}

void
InitialPlacer::ignoreSparse(int min_cluster_size)
{
     // clusterID          CellPtr
  std::map<int, std::vector<GPCell*>> clusterIDtoCells;
  std::map<int, int> macro_num_map;

  // Initialization
  for(int i = 0; i < num_cluster_; i++)
    macro_num_map[i] = 0;

  for(auto& cell : db_->movableCells() )
  {
    if(!cell->isFiller() && cell->pins().size() != 0)
    {
      clusterIDtoCells[cell->clusterID()].push_back(cell);
      if(cell->isMacro())
        macro_num_map[cell->clusterID()] += 1;
    }
  }

  // ClutserID before removal -> ClusterID after removal
  std::set<int> sparseClusterList;

  int effClusterID = 0;
  for(const auto& [cluster_id, cells] : clusterIDtoCells)
  {
    int num_cell  = static_cast<int>(cells.size());
    int num_macro = macro_num_map[cluster_id];

    // Detect Sparse Cluster
    if(num_macro == 0 && num_cell < min_cluster_size)
      sparseClusterList.insert(cluster_id);
    else
    {
      before2after_[cluster_id] = effClusterID;
      after2before_[effClusterID] = cluster_id;
      effClusterID++;
    }
  }

  garbageClusterID_ = effClusterID;
  num_cluster_      = effClusterID + 1;

  for(auto& cell : db_->movableCells() )
  {
    if(!cell->isFiller() && cell->pins().size() != 0)
    {
      int cluster_id = cell->clusterID();
      if(sparseClusterList.count(cluster_id))
        cell->setClusterID(garbageClusterID_);
      else
        cell->setClusterID(before2after_[cluster_id]);
    }
  }

  printf("Removed %d sparse clusters...\n", sparseClusterList.size());
  printf("NumCluster after Removal: %d\n", num_cluster_);

  assert(effClusterID == after2before_.size());
}

void
InitialPlacer::doSdpInit()
{
  auto build_start = getChronoNow();
  buildHostStructure();
  const double build_time = evalTime(build_start);
  printf("buildStructure       finished (takes %5.2f s)\n", build_time);

  auto cluster_start = getChronoNow();
  doClustering(host_structure_);
  cudaFreeHost(host_structure_.vertexCommunity);
  cudaFreeHost(host_structure_.communityWeight);
  cudaFreeHost(host_structure_.edges);
  cudaFreeHost(host_structure_.weights);
  cudaFreeHost(host_structure_.edgesIndex);
  cudaFreeHost(host_structure_.originalToCommunity);
  const double cluster_time = evalTime(cluster_start);
  printf("doClustering         finished (takes %5.2f s)\n", cluster_time);

  auto refine_start = getChronoNow();
  refineCluster();
  db_->setNumCluster(num_cluster_);
  computeSizeMap();

  const double refine_time = evalTime(refine_start);
  printf("refineCluster        finished (takes %5.2f s)\n", refine_time);

  auto laplacian_start = getChronoNow();

  computeFixedInfo();

  createClusterLaplacian(L_);

  extractPartialLaplacian(xf_, yf_, L_, Lff_, Lmm_, Lmf_xf_, Lmf_yf_);

  const double laplacian_time = evalTime(laplacian_start);
  printf("createLaplacian      finished (takes %5.2f s)\n", laplacian_time);

  auto sdp_start = getChronoNow();

  const double die_lx = DbToSdpX(db_->die()->lx());
  const double die_ly = DbToSdpY(db_->die()->ly());
  const double die_ux = DbToSdpX(db_->die()->ux());
  const double die_uy = DbToSdpY(db_->die()->uy());

  const double die_w = die_ux - die_lx;
  const double die_h = die_uy - die_ly;

  const double Kx = ip_param_.alpha_x * die_w * die_w * num_cluster_;
  const double Ky = ip_param_.alpha_y * die_h * die_h * num_cluster_;

  xm_ = solveSDP(Lmm_, Lmf_xf_, v_vector_, Kx);
  ym_ = solveSDP(Lmm_, Lmf_yf_, v_vector_, Ky);

  const double sdp_time = evalTime(sdp_start);
  printf("solveSDP             finished (takes %5.2f s)\n", sdp_time);

  auto export_start = getChronoNow();

  double die_lxLL = db_->die()->lx();
  double die_lyLL = db_->die()->ly();

  double die_uxLL = db_->die()->ux();
  double die_uyLL = db_->die()->uy();

  double die_cxLL = db_->die()->cx();
  double die_cyLL = db_->die()->cy();

  double mean_x = die_cxLL;  
  double mean_y = die_cyLL;  

  double deviation_x = db_->die()->dx() * ip_param_.filler_dev_coeff_x;
  double deviation_y = db_->die()->dy() * ip_param_.filler_dev_coeff_y;

  std::default_random_engine gen;
  std::normal_distribution<double> filler_noise_x(mean_x, deviation_x);
  std::normal_distribution<double> filler_noise_y(mean_y, deviation_y);

  int numMacro = db_->numMacro();

  std::normal_distribution<double> noise_to_avoid_stick_x(0.0, 200.0);
  std::normal_distribution<double> noise_to_avoid_stick_y(0.0, 200.0);

  auto& movable_cells = db_->movableCells();
  const int num_movable = movable_cells.size();

  #pragma omp parallel for num_threads(ip_param_.num_cpu_threads)
  for(int i = 0; i < num_movable; i++)
  {
    GPCell* cell = movable_cells[i];
    if(cell->isFiller())
    {
      double loc_x = filler_noise_x(gen);
      double loc_y = filler_noise_y(gen);

      loc_x = getMirrorX(loc_x, die_cxLL, die_lxLL, die_uxLL);
      loc_y = getMirrorY(loc_y, die_cyLL, die_lyLL, die_uyLL);

      cell->setCenterLocation(loc_x, loc_y);
      db_->moveCellInsideLayout(cell);
      continue;
    }

    int cID = cell->clusterID();

    double newX = SdpToDbX(xm_(cID));
    double newY = SdpToDbY(ym_(cID));

    // We have to put some noise to prevent 
    // macros in the same cluster stick together
    // ( for macro-heavy designs (e.g. bigblue2) )
    if(cell->isMacro() == true && numMacro > 500)
    {
      newX += noise_to_avoid_stick_x(gen);
      newY += noise_to_avoid_stick_y(gen);
    }

    cell->setCenterLocation( newX, newY );
    db_->moveCellInsideLayout(cell);
  }

  const double export_time = evalTime(export_start);
  printf("exportSDP            finished (takes %5.2f s)\n", export_time);
}

double
InitialPlacer::DbToSdpX(double llx) const
{
  return llx - db_->die()->cx();
}

double
InitialPlacer::DbToSdpY(double lly) const
{
  return lly - db_->die()->cy();
}

double
InitialPlacer::SdpToDbX(double ccx) const
{
  return ccx + db_->die()->cx();
}

double
InitialPlacer::SdpToDbY(double ccy) const
{
  return ccy + db_->die()->cy();
}

double
InitialPlacer::getMirrorX(double loc_x, double die_cx, double die_lx, double die_ux) const
{
  if(loc_x <= die_cx) 
    return die_lx + (die_cx - loc_x);
  else
    return die_ux - (loc_x - die_cx);
}

double
InitialPlacer::getMirrorY(double loc_y, double die_cy, double die_ly, double die_uy) const
{
  if(loc_y <= die_cy) 
    return die_ly + (die_cy - loc_y);
  else
    return die_uy - (loc_y - die_cy);
}

}; // namespace skyplace 
