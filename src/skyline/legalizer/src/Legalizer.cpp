#include <cstdio>
#include <cassert>
#include <algorithm>
#include <optional>
#include <iostream>

#include "util/Chrono.h"

#include "db/dbDatabase.h"
#include "db/dbTech.h"
#include "db/dbDesign.h"
#include "db/dbNet.h"
#include "db/dbBTerm.h"
#include "db/dbITerm.h"

#include "legalizer/Legalizer.h"

#include "Grid.h"
#include "objects/LGCell.h"
#include "objects/LGRegion.h"
#include "objects/LGVioStat.h"
#include "LGHyperParameters.h"

#include "database/CudaDatabase.h"
#include "admm_legalizer/CudaAdmmLegalizer.h"
#include "post_process/CudaPostProcess.h"

namespace legalizer 
{

Legalizer::Legalizer() {}
Legalizer::~Legalizer() {}

Legalizer::Legalizer(std::shared_ptr<db::dbDatabase> db) 
  : dbDatabase_           (db), 
    original_hpwl_        (0), 
    size_file_            (std::string()),
    constraint_           (std::string()),
    max_util_constraint_  (0.0),
    max_disp_constraint_  (0.0),
    sum_cell_area_        (0),
    num_total_rows_       (0)
{
  initParameters();
}

void
Legalizer::initParameters()
{
  params_ = std::make_shared<LGHyperParameters>();
  params_->flag_iccad17 = false;
  params_->standard_admm = false;
  params_->log_freq = 100;
  params_->num_cpu_threads = 16;
  params_->max_disp_in_row_height = 3;
  params_->max_admm_iter = 25000;
  params_->admm_terminate_thr = 0;
  params_->x_hint = 250;
  params_->y_hint = 12; //6;
  params_->iter_update_partition = 1;
  params_->max_disp_coeff = 40.0f;
  params_->tech_penalty = 1.0f;
  params_->edge_spacing_penalty = 10000.0f;

  // LEGALM style ADMM
  params_->iter_threshold_dual_update = 300;
  params_->iter_update_dual_step = 100;
  params_->iter_perturbation_on = 2;
  params_->coeff_to_init_rho = 3.0f;         // alpha_sigma (in paper)
  params_->coeff_to_init_dual = 0.5f;        // alpha_lambda (in paper)
  params_->coeff_to_update_dual_step = 0.2f; // alpha_h_f (in paper)

  // QP ADMM
  params_->qp_admm_on = false;
  params_->qp_admm_dp_site_align = true;
  params_->qp_admm_max_iter = 10000;
  params_->qp_admm_min_iter = 5;
  params_->qp_admm_rho = 2.0f;

  // Adaptive Search
  params_->use_adaptive_search = false;
  params_->stagnation_window = 10;
  params_->stagnation_threshold = 0.01f;
  params_->min_iter_to_call_adaptive_search = 300;
  params_->cooltime_to_call_adaptive_search = 100;
  params_->ovf_ratio_to_call_adaptive_search = 1.05f;
  params_->coeff_to_neglect_displace = 20.0f;

  // Sorting Policy
  params_->sorting_policy = LGSortingPolicy::OVERFLOW;
  params_->sorting_algorithm = LGSortingAlgorithm::M_BITONIC;
  params_->random_sorting_seed = 0;
}

void
Legalizer::setFlagICCAD2017(bool val)
{
  params_->flag_iccad17 = val;
}

void
Legalizer::setTechPenalty(float val)
{
  params_->tech_penalty = val;
}

void
Legalizer::setStandardAdmm(bool val)
{
  params_->standard_admm = val;
}

void
Legalizer::setQpRefine(bool val)
{
  params_->qp_admm_on = val;
}

void 
Legalizer::setCoeffToInitRho(float val)
{
  params_->coeff_to_init_rho = val;
}

void 
Legalizer::setCoeffToInitDual(float val)
{
  params_->coeff_to_init_dual = val;
}

void
Legalizer::setMaxAdmmIter(int iter)
{
  params_->max_admm_iter = iter;
}

void
Legalizer::setAdaptiveSearch(bool val)
{
  params_->use_adaptive_search = val;
}

void
Legalizer::setMinIterToCallAdaptiveSearch(int iter)
{
  params_->min_iter_to_call_adaptive_search = iter;
}

void
Legalizer::setOvfRatioToCallAdaptiveSearch(float ratio)
{
  params_->ovf_ratio_to_call_adaptive_search = ratio;
}

void
Legalizer::setCoeffToNeglectDisplace(float val)
{
  params_->coeff_to_neglect_displace = val;
}

void
Legalizer::setSortingPolicy(std::string_view policy)
{
  if(policy == "RANDOM")
    params_->sorting_policy = LGSortingPolicy::RANDOM;
  else if(policy == "SIZE")
    params_->sorting_policy = LGSortingPolicy::SIZE;
  else if(policy == "OVERFLOW")
    params_->sorting_policy = LGSortingPolicy::OVERFLOW;
  else // Default
    params_->sorting_policy = LGSortingPolicy::OVERFLOW;
}

void
Legalizer::setSortingAlgorithm(std::string_view algo)
{
  if(algo == "S_BITONIC")
    params_->sorting_algorithm = LGSortingAlgorithm::S_BITONIC;
  else if(algo == "M_BITONIC")
    params_->sorting_algorithm = LGSortingAlgorithm::M_BITONIC;
  else if(algo == "S_MERGE")
    params_->sorting_algorithm = LGSortingAlgorithm::S_MERGE;
  else if(algo == "M_MERGE")
    params_->sorting_algorithm = LGSortingAlgorithm::M_MERGE;
  else if(algo == "CUB")
    params_->sorting_algorithm = LGSortingAlgorithm::CUB;
  else // Default
    params_->sorting_algorithm = LGSortingAlgorithm::M_BITONIC;
}

void
Legalizer::setRandomSortingSeed(int seed)
{
  params_->random_sorting_seed = seed;
}

void
Legalizer::run()
{
  auto legal_start = util::getChronoNow();

  importDB();

  adjustYHint();

  makeDirectionVectors();

  detectIrregularShape();

  initialLegalize();

  initializeCudaDatabase();

  runAdmmLegalize();

  runPostProcess();

  cuda_database_->copyToHostCells(cells_);

  exportDB();

  double legal_time = util::evalTime(legal_start);

  auto lg_stat = computeStat();
  auto vio_stat = cuda_database_->checkPlace();

  printStats(lg_stat, vio_stat, legal_time);
}

void
Legalizer::adjustYHint()
{
  int max_height = 0;
  for(const auto [height, num_insts] : height_distribution_)
    max_height = std::max(max_height, height);

  int max_height_in_grid = max_height / row_height_;

  if(params_->y_hint < 3 * max_height_in_grid)
  {
    int new_y_hint = 3 * max_height_in_grid;
    printf("y_hint is adjusted to %d -> %d\n", params_->y_hint, new_y_hint);
    params_->y_hint = new_y_hint;
    // Too small y_hint can make some cells will not belong to
    // any of partitions during triplefold partitioning.
  }
}

void
Legalizer::makeDirectionVectors()
{
  // Default search space
  direction_default_.clear();
  constexpr int y_range_default = 6;
  for(int y = -y_range_default; y <= +y_range_default; y++)
  {
    int x_range = 0;
    switch(std::abs(y))
    {
      case 0: x_range = 40; break;
      case 1: x_range = 20; break;
      case 2: x_range = 10; break;
      case 3: x_range = 5;  break;
      case 4: x_range = 2;  break;
      case 5: x_range = 1;  break;
      case 6: x_range = 0;  break;
      default: break;
    }

    for(int x = -x_range; x <= +x_range; x++)
    {
      direction_default_.push_back(x);
      direction_default_.push_back(y);
    }
  }

  // Horizontal search space
  direction_horizontal_.clear();
  constexpr int y_range_horizontal = 2;
  for(int y = -y_range_horizontal; y <= +y_range_horizontal; y++)
  {
    int x_range = 0;
    switch(std::abs(y))
    {
      case 0: x_range = 24; break;
      case 1: x_range = 24; break;
      case 2: x_range = 24; break;
      default: break;
    }

    for(int x = -x_range; x <= +x_range; x++)
    {
      direction_horizontal_.push_back(x);
      direction_horizontal_.push_back(y);
    }
  }

  assert(direction_horizontal_.size() == direction_default_.size());

  // Vertical search space
  direction_vertical_.clear();
  constexpr int x_range_vertical = 2;
  for(int x = -x_range_vertical; x <= +x_range_vertical; x++)
  {
    int y_range = 0;
    switch(std::abs(x))
    {
      case 0: y_range = 24; break;
      case 1: y_range = 24; break;
      case 2: y_range = 24; break;
      default: break;
    }

    for(int y = -y_range; y <= +y_range; y++)
    {
      direction_vertical_.push_back(x);
      direction_vertical_.push_back(y);
    }
  }

  assert(direction_vertical_.size() == direction_default_.size());

  direction_wide_.clear();
  constexpr int y_range_wide = 12;
  for(int y = -y_range_wide; y <= +y_range_wide; y++)
  {
    int x_range = 0;
    switch(std::abs(y))
    {
      case  0: x_range = 280; break;
      case  1: x_range = 240; break;
      case  2: x_range = 200; break;
      case  3: x_range = 160; break;
      case  4: x_range = 120; break;
      case  5: x_range =  80; break;
      case  6: x_range =  40; break;
      case  7: x_range =  20; break;
      case  8: x_range =  10; break;
      case  9: x_range =   4; break;
      case 10: x_range =   2; break;
      case 11: x_range =   1; break;
      case 12: x_range =   0; break;
      default: break;
    }

    for(int x = -x_range; x <= +x_range; x++)
    {
      direction_wide_.push_back(x);
      direction_wide_.push_back(y);
    }
  }

  //printf("direction_wide_ %ld\n", direction_wide_.size());
}

void
Legalizer::initializeCudaDatabase()
{
  auto init_cuda = util::getChronoNow();

  // Make CudaDatabase
  cuda_database_ = std::make_shared<CudaDatabase>(regions_, params_, grid_);
  cuda_database_->copyFromHostCells(cells_, cell_tech_infos_);
  cuda_database_->copyEdgeSpacingRules(edge_spacing_rules_);

  printf("initCudaDatabase   finished (takes %5.2f s)\n", util::evalTime(init_cuda));
}

void
Legalizer::runAdmmLegalize()
{
  auto admm_start = util::getChronoNow();

  // Make CudaAdmmSolver 
  admm_legalizer::CudaAdmmLegalizer cuda_admm_legalizer(
    direction_default_, 
    direction_horizontal_, 
    direction_vertical_,
    direction_wide_,
    height_distribution_, 
    grid_, 
    params_, 
    cuda_database_);

  cuda_admm_legalizer.solve();

  printf("admmLegalize       finished (takes %5.2f s)\n", util::evalTime(admm_start));
}

void
Legalizer::runPostProcess()
{
  auto post_start = util::getChronoNow();

  CudaPostProcess cuda_post_process(params_, cuda_database_);
  cuda_post_process.doPostProcess();

  printf("postProcess        finished (takes %5.2f s)\n", util::evalTime(post_start));
}

int64_t
Legalizer::computeHpwlFromDB() const
{
  int64_t sum_hpwl = 0;

  const std::vector<dbNet*>& db_nets = dbDatabase_->getDesign()->getNets();
  for(const auto& dbnet : db_nets)
  {
    int lx = std::numeric_limits<int>::max();
    int ly = std::numeric_limits<int>::max();
    int ux = std::numeric_limits<int>::min();
    int uy = std::numeric_limits<int>::min();

    for(const auto& iterm : dbnet->getITerms())
    {
      int avg_x, avg_y;
      iterm->getAvgXY(avg_x, avg_y);
      lx = std::min(lx, avg_x);
      ly = std::min(ly, avg_y);
      ux = std::max(ux, avg_x);
      uy = std::max(uy, avg_y);
    }

    for(const auto& bterm : dbnet->getBTerms())
    {
      lx = std::min(lx, bterm->cx());
      ly = std::min(ly, bterm->cy());
      ux = std::max(ux, bterm->cx());
      uy = std::max(uy, bterm->cy());
    }

    int hpwl_this_net = (ux - lx) + (uy - ly);
    sum_hpwl += static_cast<int64_t>(hpwl_this_net);
  }

  return sum_hpwl;
}

LGStat
Legalizer::computeStat() const
{
  const auto design = dbDatabase_->getDesign();
  const double num_cell = static_cast<double>(design->getInsts().size());
  const double dbu = static_cast<double>(dbDatabase_->getTech()->getDbu());
  const std::vector<dbNet*>& db_nets = dbDatabase_->getDesign()->getNets();

  int64_t max_disp = 0.0;
  int64_t total_disp = 0.0;
  for(const auto& lg_cell : cells_)
  {
    int init_lx = lg_cell.getInitLx();
    int init_ly = lg_cell.getInitLy();

    int lg_lx = lg_cell.getLx();
    int lg_ly = lg_cell.getLy();

    int disp_x = std::abs(init_lx - lg_lx);
    int disp_y = std::abs(init_ly - lg_ly);
    int64_t disp = static_cast<int64_t>(disp_x + disp_y);

    total_disp += disp;
    max_disp = disp > max_disp ? disp : max_disp;
  }

  int64_t lg_hpwl = computeHpwlFromDB();
  
  double max_disp_double       = static_cast<double>(max_disp) / dbu;
  double total_disp_double     = static_cast<double>(total_disp) / dbu;
  double avg_disp_double       = static_cast<double>(total_disp_double) / num_cell;

  double original_hpwl_double  = static_cast<double>(original_hpwl_) / dbu;
  double lg_hpwl_double        = static_cast<double>(lg_hpwl) / dbu;

  return {max_disp_double, avg_disp_double, total_disp_double, 
          original_hpwl_double, lg_hpwl_double};
}

void
Legalizer::detectIrregularShape()
{
  for(auto region : regions_)
  {
    const int lx = region->getLx();
    const int ly = region->getLy();
    const int ux = region->getUx();
    const int uy = region->getUy();

    const int region_w = region->getWidth();
    const int region_h = region->getHeight();

    const float aspect_ratio = float(region_w) / float(region_h);

    if(aspect_ratio > 5.0f or aspect_ratio < 0.2f)
    {
      //printf("Irregular id: %d (%d, %d) - (%d, %d) aspect_ratio: %f\n", 
      //  region->getIndex(), lx, ly, ux, uy, aspect_ratio);
      if(region_w > region_h)
        region->setShape(RegionShape::THIN_HORIZONTAL);
      else
        region->setShape(RegionShape::THIN_VERTICAL);
    }
  }
}



} // namespace legalizer 
