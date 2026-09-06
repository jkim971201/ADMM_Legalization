#include "util/Chrono.h"
#include "database/CudaDatabase.h"
#include "CudaAdmmLegalizer.h"
#include "Grid.h"
#include "AdmmKernels.hpp"
#include "SortingKernels.hpp"

#include <fstream>
#include <algorithm>
#include <random>

#include <cub/cub.cuh>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/for_each.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/functional.h>

#include "cuda_linalg/CudaVectorAlgebra.h"

#include "LGHyperParameters.h"

#include "DebugHelper.hpp"
#include "GlobalUtil.h"

namespace admm_legalizer
{

void streamCompaction(
  int length,
  int target_integer,
  const cuda_linalg::CudaVector<int>& stencil, // input_data
        cuda_linalg::CudaVector<int>& output,
        cuda_linalg::CudaVector<int>& workspace)
{
  // device lambda cannot be defined within private function
  thrust::device_vector<int>::iterator iter_end =
    thrust::copy_if(
      thrust::make_counting_iterator<int>(0), thrust::make_counting_iterator<int>(length),
      stencil.begin(),                             // Stencil
      workspace.begin(),                           // Result (workspace)
      thrust::placeholders::_1 == target_integer); // Predicate
  // mysterius bug when using lambda predicator

  size_t num_target_level = thrust::distance(workspace.begin(), iter_end);

  // Copy to output
  output.resize(num_target_level);
  thrust::copy(workspace.begin(), workspace.begin() + num_target_level, output.begin());
}


CudaAdmmLegalizer::CudaAdmmLegalizer(
  const std::vector<int>& direction_default,
  const std::vector<int>& direction_horizontal,
  const std::vector<int>& direction_vertical,
  const std::vector<int>& direction_wide,
  const std::map<int, int>& height_distribution,
  const std::shared_ptr<Grid> grid,
  const std::shared_ptr<LGHyperParameters> param,
        std::shared_ptr<CudaDatabase> database) : database_(database), params_(param)
{
  num_cells_ = database_->getNumCells();

  // Copy grid information
  core_lx_ = grid->getLx();
  core_ly_ = grid->getLy();

  x_size_ = grid->getXSize();
  y_size_ = grid->getYSize();

  site_width_ = grid->getSiteWidth();
  row_height_ = grid->getRowHeight();

  pixel_dual_.resize(x_size_ * y_size_);
  pixel_color_.resize(x_size_ * y_size_);

  // Copy partition information
  num_partition_col_ = x_size_ % params_->x_hint == 0
    ? x_size_ / params_->x_hint : x_size_ / params_->x_hint + 1;
  num_partition_row_ = y_size_ % params_->y_hint == 0
    ? y_size_ / params_->y_hint : y_size_ / params_->y_hint + 1;

  // printf("NumPartCol: %d NumPartRow: %d\n", num_partition_col_, num_partition_row_);

  sequence_partition_.resize(num_partition_col_ * num_partition_row_ + 1); 
  // We will assign the last index for excluded cells
  thrust::sequence(sequence_partition_.begin(), sequence_partition_.end(), 0);

  partition_id_to_offset_.resize(sequence_partition_.size() + 1);
  partition_id_to_num_cells_.resize(sequence_partition_.size());

  // Copy list of candidate moves
  direction_default_ = direction_default;
  direction_horizontal_ = direction_horizontal;
  direction_vertical_ = direction_vertical;
  direction_wide_ = direction_wide;

  const int num_directions = direction_default.size() / 2;

  num_directions_padded_ = findMinPowerOfTwoNoLessThan(num_directions);

  // Copy height distribution
  std::vector<int> host_k_height_distribution(height_distribution.size(), 0);
  for(const auto& [height, num_cells] : height_distribution)
  {
    const int k_index = height / row_height_ - 1;
    host_k_height_distribution[k_index] = num_cells;
    assert(num_cells > 0);
  }
  k_height_to_num_cells_ = host_k_height_distribution;

  // Random Priority
  cell_random_priority_.resize(num_cells_);
  if(params_->sorting_policy == legalizer::LGSortingPolicy::RANDOM)
  {
    std::vector<int> host_random_priority(num_cells_);
    std::iota(host_random_priority.begin(), host_random_priority.end(), 0);

    std::mt19937 rand_gen(params_->random_sorting_seed);
    std::shuffle(host_random_priority.begin(), host_random_priority.end(), rand_gen);
    cell_random_priority_ = host_random_priority;
  }

  size_t shared_mem_size = getMaxSharedMemsize();
  printf("SharedMemorySize: %ld\n", shared_mem_size);

  size_t shared_mem_size_optin = getMaxSharedMemsizeOptin();
  printf("SharedMemorySizeOptin: %ld\n", shared_mem_size_optin);
}

void
CudaAdmmLegalizer::importFromGlobalDatabase()
{
  const auto& lx_grid = database_->getCellLxInGrid();
  const auto& ly_grid = database_->getCellLyInGrid();

  const size_t num_cells = num_cells_;

  cell_lx_in_grid_.resize(num_cells);
  cell_ly_in_grid_.resize(num_cells);

  thrust::copy(lx_grid.begin(), lx_grid.end(), cell_lx_in_grid_.begin());
  thrust::copy(ly_grid.begin(), ly_grid.end(), cell_ly_in_grid_.begin());

  // Resize demand-related vector
  cell_demand_plus_area_.resize(num_cells);
  cell_inverse_augmented_demand_.resize(num_cells);

  // Resize partition-related vector
  cell_id_to_partition_id_.resize(num_cells);
  cell_id_grouped_by_partition_.resize(num_cells);
  cell_key_grouped_by_partition_.resize(num_cells);

  // Resize overflow by clles vector
  cell_id_to_is_ovf_.resize(num_cells);
}

void
CudaAdmmLegalizer::updatePartition(int iter)
{
  const int num_cells = num_cells_;

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const int triplefold_number = (iter / params_->iter_update_partition) % 3;
  const int offset_x = params_->x_hint / 3 * triplefold_number;
  const int offset_y = params_->y_hint / 3 * triplefold_number;

  paintPixels(offset_x, offset_y);

  const int num_threads_per_block = 512;
  const int num_block = (num_cells + num_threads_per_block - 1) / num_threads_per_block;

  assignPartitionKernel<<<num_block, num_threads_per_block>>>(
    num_cells,
    params_->x_hint,
    params_->y_hint,
    x_size_,
    y_size_,
    offset_x,
    offset_y,
    num_partition_col_,
    num_partition_row_,
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data(),
    cell_id_to_partition_id_.data());

  /* Group cells by partition */
  // Initialize grouped indices
  thrust::sequence(cell_id_grouped_by_partition_.begin(), cell_id_grouped_by_partition_.end(), 0);

  thrust::stable_sort_by_key(
    cell_id_to_partition_id_.begin(), cell_id_to_partition_id_.end(), // Key
    cell_id_grouped_by_partition_.begin()); // Values

  // We want the last element is num_cells, so fill with num_cells when start.
  thrust::fill(partition_id_to_offset_.begin(), partition_id_to_offset_.end(), num_cells);

  thrust::lower_bound(
    cell_id_to_partition_id_.begin(), cell_id_to_partition_id_.end(), // Data
    sequence_partition_.begin(), sequence_partition_.end(), // Values to search 
    partition_id_to_offset_.begin() // Index?
  );

  // dont forget plus one
  // (this is because thrust::adjacent_difference returns x_{i} - x_{i - 1})
  thrust::adjacent_difference(
    partition_id_to_offset_.begin() + 1, partition_id_to_offset_.end(),
    partition_id_to_num_cells_.begin());

  const int num_valid_partitions = num_partition_col_ * num_partition_row_;
  max_num_cell_in_partition_ = *thrust::max_element(
    partition_id_to_num_cells_.begin(), 
    partition_id_to_num_cells_.begin() + num_valid_partitions);
  // num_cells of invalid partition can be so large that
  // it can exceed the limit of shared memory size.
}

void
CudaAdmmLegalizer::paintPixels(int offset_x, int offset_y)
{
  const int num_pixels = static_cast<int>(pixel_color_.size());
  const int num_threads_per_block = 512;
  const int num_block = (num_pixels + num_threads_per_block - 1) / num_threads_per_block;

  paintPixelsKernel<<<num_block, num_threads_per_block>>>(
    num_pixels,
    params_->x_hint,
    params_->y_hint,
    x_size_,
    y_size_,
    offset_x,
    offset_y,
    num_partition_col_,
    num_partition_row_,
    pixel_color_.data()); // color = partition_id
}

void
CudaAdmmLegalizer::solve()
{
  // Reset minimum overflow
  min_ovf_so_far_ = std::numeric_limits<int>::max();
  iter_record_.clear();
  
  // Reset last iteration when slow convergence is detected
  last_iter_slow_convergence_detected_ = -1;

  // Reset dual variables
  pixel_dual_.fillZero();

  // Fill cell_lx_in_grid_, cell_ly_in_grid_
  importFromGlobalDatabase();
  
  // Compute pixel usage (Px) 
  database_->updateUsageInGrid(cell_lx_in_grid_, cell_ly_in_grid_);
  
  // Compute pixel num_l/r_edge_type_1/2
  database_->updatePixelEdgeMapInGrid(cell_lx_in_grid_, cell_ly_in_grid_);

  // Initialize ADMM parameters
  initializeRho();
  initializeDual();

  // Main ADMM Loop
  applyAdmm();

  // Export to global database (CudaDatabase)
  exportToGlobalDatabase();
}

void
CudaAdmmLegalizer::applyAdmm()
{
  // Main ADMM Loop
  bool run_exceptional_primal_iter = false;
  for(int iter = 0; iter < params_->max_admm_iter; iter++)
  {
    if(iter % params_->iter_update_partition == 0)
      updatePartition(iter);

    if(iter % params_->iter_update_dual_step == 0 and iter > params_->iter_threshold_dual_update)
      updateDualStepSize();

    // Step 1: Primal Update
    if(run_exceptional_primal_iter == false or params_->use_adaptive_search == false)
      doPrimalUpdate(iter);
    else
      doPrimalUpdateForSlowConvergence(iter);

    // Step 2: Dual Update
    doDualUpdate();

    const auto& [ovf, disp] = examineIteration();

    run_exceptional_primal_iter = detectSlowConvergence(iter, ovf, min_ovf_so_far_);

    bool convergence = ovf <= params_->admm_terminate_thr ? true : false;
    if(convergence == true or iter % params_->log_freq == 0)
      printProgress(iter, ovf, disp);

    if(convergence == true)
    {
      doFinalRefine();
      break;
    }
  }
}

void
CudaAdmmLegalizer::initializeRho()
{
  auto& pixel_usage = database_->getPixelUsage();

  const int num_cells = num_cells_;
  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  int num_threads_per_block = 512;
  int num_block = (num_cells + num_threads_per_block - 1) / num_threads_per_block;

  computeDemandKernel<<<num_block, num_threads_per_block>>>(
    num_cells,
    x_size_,
    y_size_,
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    pixel_usage.data(),
    cell_demand_plus_area_.data());

  float demand_sum = computeVectorSum(cell_demand_plus_area_);

  float initial_rho = params_->coeff_to_init_rho * demand_sum / float(num_cells);

  rho_ = initial_rho;
  printf("Initial Rho : %f (DemandSum: %f)\n", initial_rho, demand_sum);
}

float
CudaAdmmLegalizer::computeInverseAugmentedDemand()
{
  const auto& pixel_usage = database_->getPixelUsage();

  const int num_cells = num_cells_;
  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const float expected_avg_disp = float(row_height_) / float(site_width_);

  int num_threads_per_block = 512;
  int num_block = (num_cells + num_threads_per_block - 1) / num_threads_per_block;

  computeInverseAugmentedDemandKernel<<<num_block, num_threads_per_block>>>(
    num_cells,
    x_size_,
    y_size_,
    rho_,
    expected_avg_disp,
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    pixel_usage.data(),
    cell_inverse_augmented_demand_.data());

   float sum_inverse_augmented_demand = computeVectorSum(cell_inverse_augmented_demand_);
   return sum_inverse_augmented_demand;
}

void
CudaAdmmLegalizer::initializeDual()
{
  float sum_inverse_augmented_demand = computeInverseAugmentedDemand();

  float initial_dual = params_->coeff_to_init_dual * sum_inverse_augmented_demand / float(num_cells_);

  thrust::fill(pixel_dual_.begin(), pixel_dual_.end(), initial_dual);
  initial_dual_ = initial_dual;
  dual_step_size_ = initial_dual;

  printf("Initial Dual: %f (Inverse Augmented Sum: %f)\n", initial_dual, sum_inverse_augmented_demand);
}

void
CudaAdmmLegalizer::computeKeyArrayForSorting(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
  const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
  const cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition,
        cuda_linalg::CudaVector<int>& cell_key_grouped_by_partition)
{
  const int num_threads_per_block = 256;
  const int num_block = (num_cells_ + num_threads_per_block - 1) / num_threads_per_block;

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();
  const auto& pixel_usage = database_->getPixelUsage();

  int coeff_random = 0;
  int coeff_ovf = 10;
  int coeff_area = 1;
  if(params_->sorting_policy == legalizer::LGSortingPolicy::RANDOM)
  {
    coeff_random = 1;
    coeff_ovf = 0;
    coeff_area = 0;
  }

  computeKeyArrayForSortingKernel<<<num_block, num_threads_per_block>>>(
    x_size_,
    y_size_,
    num_cells_,
    coeff_random,
    coeff_ovf,
    coeff_area,
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    lx_in_grid.data(),
    ly_in_grid.data(),
    pixel_usage.data(),
    partition_id_to_num_cell.data(),
    partition_id_to_offset.data(),
    cell_random_priority_.data(),
    cell_id_grouped_by_partition.data(),
    cell_key_grouped_by_partition.data());
}

void
CudaAdmmLegalizer::doBitonicSort(
  const int num_threads_per_block,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_ovf_grouped_by_partition, // Key
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)  // Value
{
  const int array_size_padded = findMinPowerOfTwoNoLessThan(max_num_cell_in_partition_);
  const size_t shared_mem_size = 2 * sizeof(int) * array_size_padded;
  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  bitonicSortKernel<<<num_valid_partitions, num_threads_per_block, shared_mem_size>>>(
    array_size_padded, /* buffer size */
    partition_id_to_offset.data(),
    cell_id_grouped_by_partition.data(),
    cell_key_grouped_by_partition_.data());
}

void
CudaAdmmLegalizer::doMergeSort(
  const int num_threads_per_block,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_ovf_grouped_by_partition, // Key
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)  // Value
{
  const int array_size_padded = findMinPowerOfTwoNoLessThan(max_num_cell_in_partition_);
  const size_t shared_mem_size = 4 * sizeof(int) * array_size_padded;
  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  mergeSortKernel<<<num_valid_partitions, num_threads_per_block, shared_mem_size>>>(
    array_size_padded, /* buffer size */
    partition_id_to_offset.data(),
    cell_id_grouped_by_partition.data(),
    cell_key_grouped_by_partition_.data());
}

void
CudaAdmmLegalizer::doCubSort(
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_key_grouped_by_partition, // Key
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)  // Value
{
  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  void* temp_storage = nullptr;
  size_t temp_storage_bytes = 0;

  const int* offsets = partition_id_to_offset.data();

  int* keys_in = cell_key_grouped_by_partition.data();
  int* keys_out = cell_key_grouped_by_partition.data();

  int* vals_in = cell_id_grouped_by_partition.data();
  int* vals_out = cell_id_grouped_by_partition.data();

  // Compute temp storage size
  cub::DeviceSegmentedRadixSort::SortPairsDescending(
    temp_storage,
    temp_storage_bytes,
    keys_in, 
    keys_out, 
    vals_in, 
    vals_out,
    num_cells_, 
    num_valid_partitions, 
    offsets, 
    offsets + 1);

  cudaMalloc(&temp_storage, temp_storage_bytes);

  // Real Sorting
  cub::DeviceSegmentedRadixSort::SortPairsDescending(
    temp_storage,
    temp_storage_bytes,
    keys_in, 
    keys_out, 
    vals_in, 
    vals_out,
    num_cells_, 
    num_valid_partitions,
    offsets, 
    offsets + 1);

  cudaFree(temp_storage);
}

void
CudaAdmmLegalizer::executeSorting(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
  const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)
{
  computeKeyArrayForSorting(
    lx_in_grid,
    ly_in_grid,
    partition_id_to_num_cell,
    partition_id_to_offset,
    cell_id_grouped_by_partition,
    cell_key_grouped_by_partition_);

  if(params_->sorting_algorithm == legalizer::LGSortingAlgorithm::S_MERGE)
    doMergeSort(1, partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
  else if(params_->sorting_algorithm == legalizer::LGSortingAlgorithm::M_MERGE)
    doMergeSort(32, partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
  else if(params_->sorting_algorithm == legalizer::LGSortingAlgorithm::S_BITONIC)
    doBitonicSort(1, partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
  else if(params_->sorting_algorithm == legalizer::LGSortingAlgorithm::M_BITONIC)
    doBitonicSort(32, partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
  else if(params_->sorting_algorithm == legalizer::LGSortingAlgorithm::CUB)
    doCubSort(partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
  else
    doBitonicSort(32, partition_id_to_offset, cell_key_grouped_by_partition_, cell_id_grouped_by_partition);
}

void
CudaAdmmLegalizer::doPrimalUpdate(int iter)
{
  // Sort by sum of overflow for each cell in the partition
  executeSorting(
    cell_lx_in_grid_,
    cell_ly_in_grid_,
    partition_id_to_num_cells_,
    partition_id_to_offset_,
    cell_id_grouped_by_partition_);

  const int num_cells = num_cells_;
  const int num_directions = static_cast<int>(direction_default_.size() / 2);

  const auto& pin_id_to_bgn = database_->getPinBgn();
  const auto& pin_id_to_end = database_->getPinEnd();
  const auto& pin_id_to_layer = database_->getPinLayer();
  const auto& macro_id_to_pin_offsets = database_->getPinOffsets();

  const auto& cell_macro_id = database_->getCellMacroID();

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const auto& cell_original_lx_in_dbu = database_->getCellOriginalLxInDbu();
  const auto& cell_original_ly_in_dbu = database_->getCellOriginalLyInDbu();
  const auto& cell_group_id = database_->getCellGroupId();
  const auto& cell_has_ground_at_bottom = database_->getCellHasGroundAtBottom();

  const auto& cell_id_to_ledge_type = database_->getCellLEdgeType();
  const auto& cell_id_to_redge_type = database_->getCellREdgeType();
  const auto& edge_spacing_rules = database_->getEdgeSpacingRules();

  const auto& grid_y_to_is_vdd_up = database_->getGridYToIsVddUp();

  const auto& pixel_power_metal_layer = database_->getPixelPowerMetalLayer();

  const auto& pixel_valid     = database_->getPixelValid();
  const auto& pixel_group_id  = database_->getPixelGroupId();
  const auto& pixel_region_id = database_->getPixelRegionId();
  const auto& pixel_shape     = database_->getPixelShape();
        auto& pixel_usage     = database_->getPixelUsage();

        auto& pixel_num_ledge_type1 = database_->getPixelNumLEdgeType1();
        auto& pixel_num_ledge_type2 = database_->getPixelNumLEdgeType2();
        auto& pixel_num_redge_type1 = database_->getPixelNumREdgeType1();
        auto& pixel_num_redge_type2 = database_->getPixelNumREdgeType2();

  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  const int is_standard_admm = params_->standard_admm == true ? 1 : 0;

  const int num_block = num_valid_partitions;
  const int num_threads_per_block = 256;
  const size_t shared_mem_size = 2 * sizeof(float) * num_directions_padded_;

  //const bool perturbation_on = 
  //  (iter < params_->iter_perturbation_on or iter == params_->iter_threshold_dual_update) ? true : false;

  const bool perturbation_on = (iter < params_->iter_perturbation_on) ? true : false;

  assert(num_directions_padded_ <= num_threads_per_block);

  doPrimalUpdateKernel<<<num_block, num_threads_per_block, shared_mem_size>>>(
    params_->flag_iccad17,
    is_standard_admm,
    perturbation_on,
    iter,
    num_cells,
    num_directions,
    num_directions_padded_,
    x_size_,
    y_size_,
    core_lx_,
    core_ly_,
    site_width_,
    row_height_,
    params_->max_disp_in_row_height,
    params_->max_disp_coeff,
    params_->tech_penalty,
    params_->edge_spacing_penalty,
    rho_,
    partition_id_to_offset_.data(),
    direction_default_.data(),
    direction_horizontal_.data(),
    direction_vertical_.data(),
    k_height_to_num_cells_.data(),
    pin_id_to_bgn.data(),
    pin_id_to_end.data(),
    pin_id_to_layer.data(),
    macro_id_to_pin_offsets.data(),
    cell_macro_id.data(),
    cell_id_grouped_by_partition_.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    cell_original_lx_in_dbu.data(),
    cell_original_ly_in_dbu.data(),
    cell_group_id.data(),
    cell_has_ground_at_bottom.data(),
    cell_id_to_ledge_type.data(),
    cell_id_to_redge_type.data(),
    edge_spacing_rules.data(),
    grid_y_to_is_vdd_up.data(),
    pixel_color_.data(),
    pixel_valid.data(),
    pixel_group_id.data(),
    pixel_shape.data(),
    pixel_power_metal_layer.data(),
    pixel_dual_.data(),
    pixel_usage.data(),
    pixel_num_ledge_type1.data(),
    pixel_num_ledge_type2.data(),
    pixel_num_redge_type1.data(),
    pixel_num_redge_type2.data(),
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data());
}

void
CudaAdmmLegalizer::doDualUpdate()
{
  const auto& pixel_usage = database_->getPixelUsage();

  auto zip_bgn = thrust::make_zip_iterator(thrust::make_tuple(pixel_dual_.begin(), pixel_usage.begin()));
  auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(pixel_dual_.end(), pixel_usage.end()));

  thrust::for_each(zip_bgn, zip_end, dual_update_functor_with_step_size(dual_step_size_, rho_));
}

void
CudaAdmmLegalizer::updateDualStepSize()
{
  dual_step_size_ += params_->coeff_to_update_dual_step * initial_dual_;
}

void
CudaAdmmLegalizer::sortEachPartition(
  const int num_threads_per_block,
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
  const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)
{
  const int num_cells = num_cells_;

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const auto& pixel_usage = database_->getPixelUsage();

  int array_size_padded = findMinPowerOfTwoNoLessThan(max_num_cell_in_partition_);
  assert(array_size_padded < 12288);

  const size_t shared_mem_size = 2 * sizeof(int) * array_size_padded;
  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  int coeff_random = 0;
  int coeff_ovf = 10;
  int coeff_area = 1;
  if(params_->sorting_policy == legalizer::LGSortingPolicy::RANDOM)
  {
    coeff_random = 1;
    coeff_ovf = 0;
    coeff_area = 0;
  }

  sortEachPartitionKernel<<<num_valid_partitions, num_threads_per_block, shared_mem_size>>>(
    num_cells,
    x_size_,
    y_size_,
    array_size_padded, /* buffer size */
    coeff_random,
    coeff_ovf,
    coeff_area,
    partition_id_to_offset.data(),
    pixel_usage.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    lx_in_grid.data(),
    ly_in_grid.data(),
    cell_random_priority_.data(),
    cell_id_grouped_by_partition.data());
}

void
CudaAdmmLegalizer::sortEachPartitionWithMergeSort(
  const int num_threads_per_block,
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
  const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
  const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition)
{
  const int num_cells = num_cells_;

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const auto& pixel_usage = database_->getPixelUsage();

  int array_size_padded = findMinPowerOfTwoNoLessThan(max_num_cell_in_partition_);
  assert(array_size_padded < 12288);

  const size_t shared_mem_size = 4 * sizeof(int) * array_size_padded;
  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  int coeff_random = 0;
  int coeff_ovf = 10;
  int coeff_area = 1;
  if(params_->sorting_policy == legalizer::LGSortingPolicy::RANDOM)
  {
    coeff_random = 1;
    coeff_ovf = 0;
    coeff_area = 0;
  }

  sortEachPartitionWithMergeSortKernel<<<num_valid_partitions, num_threads_per_block, shared_mem_size>>>(
    num_cells,
    x_size_,
    y_size_,
    array_size_padded, /* buffer size */
    coeff_random,
    coeff_ovf,
    coeff_area,
    partition_id_to_offset.data(),
    pixel_usage.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    lx_in_grid.data(),
    ly_in_grid.data(),
    cell_random_priority_.data(),
    cell_id_grouped_by_partition.data());
}

void
CudaAdmmLegalizer::exportToGlobalDatabase() 
{
  database_->commitPlaceInGrid(cell_lx_in_grid_, cell_ly_in_grid_);
}

bool
CudaAdmmLegalizer::detectSlowConvergence(int iter, int ovf, int min_ovf_so_far) 
{
  // Do not run exceptional primal_update too early
  if(iter < params_->min_iter_to_call_adaptive_search)
    return false;
  else if(min_ovf_so_far == 0)
    return false;
  else if(last_iter_slow_convergence_detected_ > 0 and 
          iter < last_iter_slow_convergence_detected_ + params_->cooltime_to_call_adaptive_search)
    return false; // Avoid frequent call
  else
  {
    bool detect_by_ovf_ratio = false;
    float ratio = float(ovf) / float(min_ovf_so_far);
    if(ratio >= params_->ovf_ratio_to_call_adaptive_search)
    {
      detect_by_ovf_ratio = true;
    }

    bool detect_by_stagnation = true;
    if(ovf_dequeue_.size() == params_->stagnation_window)
    {
      int ovf_front = ovf_dequeue_.front();
      for(const int ovf_dequeue : ovf_dequeue_)
      {
        if(ovf_dequeue != ovf_front)
        {
          detect_by_stagnation = false;
          break;
        }
      }
    }

    bool detect_slow_converge = detect_by_ovf_ratio or detect_by_stagnation;

    last_iter_slow_convergence_detected_ 
      = detect_slow_converge == true ? iter : last_iter_slow_convergence_detected_;

    return detect_slow_converge;
  }
}

void
CudaAdmmLegalizer::printProgress(int iter, int ovf, float disp) const
{
  const float avg_disp_in_site = disp / float(num_cells_) / float(site_width_);
  printf("Iter: %5d Ovf: %7d Disp: %10.0f (%5.2f sites) MinOvf: %3d\n", 
    iter, ovf, disp, avg_disp_in_site, min_ovf_so_far_);
}

void
CudaAdmmLegalizer::writePartition(int iter) const
{
  // Be careful that the order of cell_id_to_partition is messed up
  // after doing thrust::stable_sort_by_key
  std::vector<int> host_cell_id_to_lx(cell_lx_in_grid_.size());
  thrust::copy(cell_lx_in_grid_.begin(), cell_lx_in_grid_.end(),
    host_cell_id_to_lx.begin());

  std::vector<int> host_cell_id_to_ly(cell_ly_in_grid_.size());
  thrust::copy(cell_ly_in_grid_.begin(), cell_ly_in_grid_.end(),
    host_cell_id_to_ly.begin());

  std::vector<int> host_cell_id_to_partition_id(cell_id_to_partition_id_.size());
  thrust::copy(cell_id_to_partition_id_.begin(), cell_id_to_partition_id_.end(), 
    host_cell_id_to_partition_id.begin());

  std::string file_name = "partition_" + std::to_string(iter) + ".txt";
  std::ofstream output(file_name);

  for(int i = 0; i < host_cell_id_to_partition_id.size(); i++)
  {
    output << host_cell_id_to_lx[i] << " ";
    output << host_cell_id_to_ly[i] << " ";
    output << host_cell_id_to_partition_id[i] << " " << std::endl;
  }
}

void
CudaAdmmLegalizer::doPrimalUpdateForSlowConvergence(int iter)
{
  cell_id_to_is_ovf_.fillZero();

  cell_compaction_workspace_.resize(num_cells_);

  int num_threads_per_block = 1024;
  int num_block = (num_cells_ + num_threads_per_block - 1) / num_threads_per_block;

  auto& pixel_usage = database_->getPixelUsage();

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  extractCongestedCellsKernel<<<num_block, num_threads_per_block>>>(
    num_cells_,
    y_size_,
    pixel_usage.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data(),
    cell_id_to_is_ovf_.data());

  streamCompaction(
    num_cells_,
    1,
    cell_id_to_is_ovf_,
    ovf_cells_,
    cell_compaction_workspace_);

  //cuda_linalg::printVector(ovf_cells_, "Ovf Cell IDs");

  const auto& cell_original_lx_in_dbu = database_->getCellOriginalLxInDbu();
  const auto& cell_original_ly_in_dbu = database_->getCellOriginalLyInDbu();
  const auto& cell_group_id = database_->getCellGroupId();
  const auto& cell_has_ground_at_bottom = database_->getCellHasGroundAtBottom();

  const auto& cell_id_to_ledge_type = database_->getCellLEdgeType();
  const auto& cell_id_to_redge_type = database_->getCellREdgeType();

  const auto& grid_y_to_is_vdd_up = database_->getGridYToIsVddUp();

  const auto& pixel_valid     = database_->getPixelValid();
  const auto& pixel_group_id  = database_->getPixelGroupId();
  const auto& pixel_region_id = database_->getPixelRegionId();
  const auto& pixel_shape     = database_->getPixelShape();

        auto& pixel_num_ledge_type1 = database_->getPixelNumLEdgeType1();
        auto& pixel_num_ledge_type2 = database_->getPixelNumLEdgeType2();
        auto& pixel_num_redge_type1 = database_->getPixelNumREdgeType1();
        auto& pixel_num_redge_type2 = database_->getPixelNumREdgeType2();

  const auto& region_grid_x_min = database_->getRegionGridXMin();
  const auto& region_grid_y_min = database_->getRegionGridYMin();
  const auto& region_grid_x_max = database_->getRegionGridXMax();
  const auto& region_grid_y_max = database_->getRegionGridYMax();

  const int num_ovf_cells = static_cast<int>(ovf_cells_.size());
  const int num_directions = static_cast<int>(direction_wide_.size() / 2);
  const int num_directions_padded = findMinPowerOfTwoNoLessThan(num_directions);

  //assert(num_directions == 1024);
  assert(num_directions_padded == 4096);

  const size_t shared_mem_size = 2 * sizeof(float) * num_directions_padded;

  projectionAwareSearchKernel<<<1, num_threads_per_block, shared_mem_size>>>(
    params_->flag_iccad17,
    num_cells_,
    num_ovf_cells,
    num_directions,
    num_directions_padded,
    x_size_,
    y_size_,
    core_lx_,
    core_ly_,
    site_width_,
    row_height_,
    params_->max_disp_in_row_height,
    params_->max_disp_coeff,
    rho_,
    params_->coeff_to_neglect_displace,
    ovf_cells_.data(),
    direction_wide_.data(),
    k_height_to_num_cells_.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    cell_original_lx_in_dbu.data(),
    cell_original_ly_in_dbu.data(),
    cell_group_id.data(),
    cell_has_ground_at_bottom.data(),
    cell_id_to_ledge_type.data(),
    cell_id_to_redge_type.data(),
    grid_y_to_is_vdd_up.data(),
    region_grid_x_min.data(),
    region_grid_y_min.data(),
    region_grid_x_max.data(),
    region_grid_y_max.data(),
    pixel_color_.data(),
    pixel_valid.data(),
    pixel_group_id.data(),
    pixel_shape.data(),
    pixel_region_id.data(),
    pixel_dual_.data(),
    pixel_usage.data(),
    pixel_num_ledge_type1.data(),
    pixel_num_ledge_type2.data(),
    pixel_num_redge_type1.data(),
    pixel_num_redge_type2.data(),
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data());

  printf("Primal Update with Adaptive Search (Iter: %4d OvfCells: %3ld)\n", iter, ovf_cells_.size());
}

std::pair<int, float> 
CudaAdmmLegalizer::examineIteration()
{
  int ovf = database_->computeOverflow();
  float disp = database_->computeDisplace(cell_lx_in_grid_, cell_ly_in_grid_);
  min_ovf_so_far_ = std::min(ovf, min_ovf_so_far_);

  iter_record_.push_back({ovf, disp});

  ovf_dequeue_.push_back(ovf);
  if(ovf_dequeue_.size() > params_->stagnation_window)
    ovf_dequeue_.pop_front();

  return {ovf, disp};
}

void
CudaAdmmLegalizer::doFinalRefine()
{
  const int num_cells = num_cells_;
  const int num_directions = static_cast<int>(direction_default_.size() / 2);

  const auto& pin_id_to_bgn = database_->getPinBgn();
  const auto& pin_id_to_end = database_->getPinEnd();
  const auto& pin_id_to_layer = database_->getPinLayer();
  const auto& macro_id_to_pin_offsets = database_->getPinOffsets();

  const auto& cell_macro_id = database_->getCellMacroID();

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();

  const auto& cell_original_lx_in_dbu = database_->getCellOriginalLxInDbu();
  const auto& cell_original_ly_in_dbu = database_->getCellOriginalLyInDbu();
  const auto& cell_group_id = database_->getCellGroupId();
  const auto& cell_has_ground_at_bottom = database_->getCellHasGroundAtBottom();

  const auto& cell_id_to_ledge_type = database_->getCellLEdgeType();
  const auto& cell_id_to_redge_type = database_->getCellREdgeType();
  const auto& edge_spacing_rules = database_->getEdgeSpacingRules();

  const auto& grid_y_to_is_vdd_up = database_->getGridYToIsVddUp();

  const auto& pixel_power_metal_layer = database_->getPixelPowerMetalLayer();

  const auto& pixel_valid     = database_->getPixelValid();
  const auto& pixel_group_id  = database_->getPixelGroupId();
  const auto& pixel_region_id = database_->getPixelRegionId();
  const auto& pixel_shape     = database_->getPixelShape();
        auto& pixel_usage     = database_->getPixelUsage();

        auto& pixel_num_ledge_type1 = database_->getPixelNumLEdgeType1();
        auto& pixel_num_ledge_type2 = database_->getPixelNumLEdgeType2();
        auto& pixel_num_redge_type1 = database_->getPixelNumREdgeType1();
        auto& pixel_num_redge_type2 = database_->getPixelNumREdgeType2();

  const int num_valid_partitions = num_partition_col_ * num_partition_row_;

  const int is_standard_admm = params_->standard_admm == true ? 1 : 0;

  const int num_block = num_valid_partitions;
  const int num_threads_per_block = 512;
  const size_t shared_mem_size = 2 * sizeof(float) * num_directions_padded_;

  assert(num_directions_padded_ <= num_threads_per_block);

  doFinalRefineKernel<<<num_block, num_threads_per_block, shared_mem_size>>>(
    params_->flag_iccad17,
    num_cells,
    num_directions,
    num_directions_padded_,
    x_size_,
    y_size_,
    core_lx_,
    core_ly_,
    site_width_,
    row_height_,
    params_->max_disp_in_row_height,
    params_->max_disp_coeff,
    params_->tech_penalty,
    params_->edge_spacing_penalty,
    rho_,
    partition_id_to_offset_.data(),
    direction_default_.data(),
    direction_horizontal_.data(),
    direction_vertical_.data(),
    k_height_to_num_cells_.data(),
    pin_id_to_bgn.data(),
    pin_id_to_end.data(),
    pin_id_to_layer.data(),
    macro_id_to_pin_offsets.data(),
    cell_macro_id.data(),
    cell_id_grouped_by_partition_.data(),
    cell_w_in_grid.data(),
    cell_h_in_grid.data(),
    cell_original_lx_in_dbu.data(),
    cell_original_ly_in_dbu.data(),
    cell_group_id.data(),
    cell_has_ground_at_bottom.data(),
    cell_id_to_ledge_type.data(),
    cell_id_to_redge_type.data(),
    edge_spacing_rules.data(),
    grid_y_to_is_vdd_up.data(),
    pixel_color_.data(),
    pixel_valid.data(),
    pixel_group_id.data(),
    pixel_shape.data(),
    pixel_power_metal_layer.data(),
    pixel_dual_.data(),
    pixel_usage.data(),
    pixel_num_ledge_type1.data(),
    pixel_num_ledge_type2.data(),
    pixel_num_redge_type1.data(),
    pixel_num_redge_type2.data(),
    cell_lx_in_grid_.data(),
    cell_ly_in_grid_.data());

  //int ovf = database_->computeOverflow();
  //float disp = database_->computeDisplace(cell_lx_in_grid_, cell_ly_in_grid_);
  //printf("(disp, ovf) after final_refine: %.1f, ovf: %d\n", disp, ovf);
}

void
CudaAdmmLegalizer::writeRecord() const
{
  std::string file_name = "convergence_record.txt";
  std::ofstream output(file_name);

  for(int i = 0; i < iter_record_.size(); i++)
  {
    output << i << " ";
    output << iter_record_[i].first << " ";
    output << iter_record_[i].second << " " << std::endl;
  }
}

}
