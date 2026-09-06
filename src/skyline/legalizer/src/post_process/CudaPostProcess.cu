#include "util/Chrono.h"
#include "database/CudaDatabase.h"
#include "CudaPostProcess.h"
#include "CudaQpSolver.h"
#include "PostProcessKernel.hpp"
#include "GlobalUtil.h"
#include "LGHyperParameters.h"

#include "cuda_linalg/CudaVectorAlgebra.h"

#include <thrust/sort.h>
#include <thrust/count.h>
#include <thrust/sequence.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>

namespace legalizer
{

CudaPostProcess::CudaPostProcess(
  const std::shared_ptr<LGHyperParameters> params,
        std::shared_ptr<CudaDatabase> database)
  : database_(database), params_(params)
{}

void
CudaPostProcess::doPostProcess()
{
  // 1. Rotation
  rotateCells();

  if(params_->qp_admm_on == false)
    return;

  auto qp_refine_start = util::getChronoNow();
  float displace_before_refine = database_->computeDisplace();

  // 2. Stamp Cells
  stampCells();

  // 3. Create cut rows
  makeCutRows();

  // 4. Insert cells to cutrows
  insertCellsToCutRows();

  // 5. Sort cells by lx
  sortCellsEachCutRow();

  // 6. Solve QP
  solveQP();

  float displace_after_refine = database_->computeDisplace();

  const double qp_refine_time = util::evalTime(qp_refine_start);

  float improvement_ratio 
    = (displace_after_refine - displace_before_refine) / displace_before_refine * 100.0f;

  printf("Displacement %.1f -> %.1f (%.2f %%) Time: %.2f s\n",
    displace_before_refine,
    displace_after_refine,
    improvement_ratio,
    qp_refine_time);
}

void
CudaPostProcess::rotateCells()
{
  const auto& grid_y_to_is_vdd_up = database_->getGridYToIsVddUp();

  const auto& cell_id_to_has_ground_at_bottom = database_->getCellHasGroundAtBottom();

  const auto& cell_id_to_ly_in_grid = database_->getCellLyInGrid();
  const auto& cell_id_to_height_in_grid = database_->getCellHeightInGrid();

        auto& cell_id_to_need_rotation = database_->getCellNeedRotation();

  const int num_cells = database_->getNumCells();

  int num_thread_per_block = 512;
  int num_block = (num_cells + num_thread_per_block - 1) / num_thread_per_block;

  rotateKernel<<<num_block, num_thread_per_block>>>(
    num_cells,
    grid_y_to_is_vdd_up.data(),
    cell_id_to_ly_in_grid.data(),
    cell_id_to_height_in_grid.data(),
    cell_id_to_has_ground_at_bottom.data(),
    cell_id_to_need_rotation.data());
}

void
CudaPostProcess::stampCells()
{
  const int num_cells = database_->getNumCells();
  const int num_pixels = database_->getNumPixels();

  pixel_id_to_cell_id_.resize(num_pixels);
  thrust::fill(pixel_id_to_cell_id_.begin(), pixel_id_to_cell_id_.end(), -1);
  // -1 initialization is important!!!

  const auto& cell_id_to_lx_in_grid = database_->getCellLxInGrid();
  const auto& cell_id_to_ly_in_grid = database_->getCellLyInGrid();
  const auto& cell_id_to_width_in_grid = database_->getCellWidthInGrid();
  const auto& cell_id_to_height_in_grid = database_->getCellHeightInGrid();

  int num_thread_per_block = 512;
  int num_block = (num_cells + num_thread_per_block - 1) / num_thread_per_block;

  stampCellsKernel<<<num_block, num_thread_per_block>>>(
    num_cells,
    database_->getXSize(),
    database_->getYSize(),
    cell_id_to_lx_in_grid.data(),
    cell_id_to_ly_in_grid.data(),
    cell_id_to_width_in_grid.data(),
    cell_id_to_height_in_grid.data(),
    pixel_id_to_cell_id_.data());
}

void
CudaPostProcess::makeCutRows()
{
  // Mark Multi-height cells
  const int num_cells  = database_->getNumCells();
  const int num_pixels = database_->getNumPixels();
  const int x_size     = database_->getXSize();
  const int y_size     = database_->getYSize();

  const auto& cell_lx_in_grid = database_->getCellLxInGrid();
  const auto& cell_ly_in_grid = database_->getCellLyInGrid();
  const auto& cell_w_in_grid  = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid  = database_->getCellHeightInGrid();

  const auto& pixel_valid    = database_->getPixelValid();
  const auto& pixel_group_id = database_->getPixelGroupId();

  grid_y_to_num_cut_rows_.resize(y_size);

  int num_thread_per_block = 512;
  int num_block_grid_y = (y_size + num_thread_per_block - 1) / num_thread_per_block;

  countCutRowsKernel<<<num_block_grid_y, num_thread_per_block>>>(
    x_size, 
    y_size, 
    cell_h_in_grid.data(),
    pixel_group_id.data(),
    pixel_valid.data(),
    pixel_id_to_cell_id_.data(),
    grid_y_to_num_cut_rows_.data());

  grid_y_to_cut_row_offset_.resize(y_size + 1);

  thrust::inclusive_scan(
    grid_y_to_num_cut_rows_.begin(), grid_y_to_num_cut_rows_.end(),
    grid_y_to_cut_row_offset_.begin() + 1);

  // Get the last value of grid_y_to_cut_row_offset_
  num_cut_rows_ = 0;
  cudaMemcpy(
    &num_cut_rows_, 
    grid_y_to_cut_row_offset_.data() + grid_y_to_cut_row_offset_.size() - 1, 
    sizeof(int), 
    cudaMemcpyDeviceToHost);

  cut_row_id_to_grid_y_.resize(num_cut_rows_);
  cut_row_id_to_grid_x_bgn_.resize(num_cut_rows_);
  cut_row_id_to_grid_x_end_.resize(num_cut_rows_);
  cut_row_id_to_cell_offset_.resize(num_cut_rows_ + 2);

  generateCutRowsKernel<<<num_block_grid_y, num_thread_per_block>>>(
    x_size, 
    y_size, 
    cell_h_in_grid.data(),
    pixel_group_id.data(),
    pixel_valid.data(),
    pixel_id_to_cell_id_.data(),
    grid_y_to_cut_row_offset_.data(),
    cut_row_id_to_grid_y_.data(),
    cut_row_id_to_grid_x_bgn_.data(),
    cut_row_id_to_grid_x_end_.data());

  int num_block_cut_row 
    = (num_cut_rows_ + num_thread_per_block - 1) / num_thread_per_block;

  pixel_id_to_cut_row_id_.resize(num_pixels);
  thrust::fill(pixel_id_to_cut_row_id_.begin(), pixel_id_to_cut_row_id_.end(), -1);
  // Initialize with -1

  markCutRowToGridKernel<<<num_block_cut_row, num_thread_per_block>>>(
    x_size, 
    y_size, 
    num_cut_rows_,
    cut_row_id_to_grid_y_.data(),
    cut_row_id_to_grid_x_bgn_.data(),
    cut_row_id_to_grid_x_end_.data(),
    pixel_id_to_cut_row_id_.data());
}

void
CudaPostProcess::insertCellsToCutRows()
{
  const int num_cells = database_->getNumCells();

  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();
  const auto& cell_lx_in_grid = database_->getCellLxInGrid();
  const auto& cell_ly_in_grid = database_->getCellLyInGrid();
  
  int num_thread_per_block = 512;
  int num_block = (num_cells + num_thread_per_block - 1) / num_thread_per_block;
  cell_id_to_cut_row_id_.resize(num_cells);

  markCutRowToCellKernel<<<num_block, num_thread_per_block>>>(
    num_cells,
    num_cut_rows_, // including invalid cut_row
    database_->getXSize(),
    database_->getYSize(),
    pixel_id_to_cut_row_id_.data(),
    cell_h_in_grid.data(),
    cell_lx_in_grid.data(),
    cell_ly_in_grid.data(),
    cell_id_to_cut_row_id_.data());

  //printVector(cell_id_to_cut_row_id_, "CellIDtoCutRowID");

  cuda_linalg::CudaVector<int> partition_sequence(num_cut_rows_ + 1);
  thrust::sequence(partition_sequence.begin(), partition_sequence.end(), 0);

  // Initialize grouped indices
  cell_grouped_by_cut_row_.resize(num_cells);
  thrust::sequence(cell_grouped_by_cut_row_.begin(), cell_grouped_by_cut_row_.end(), 0);

  thrust::stable_sort_by_key(
    cell_id_to_cut_row_id_.begin(), cell_id_to_cut_row_id_.end(), // Key
    cell_grouped_by_cut_row_.begin()); // Values

  thrust::lower_bound(
    cell_id_to_cut_row_id_.begin(), cell_id_to_cut_row_id_.end(), // Data
    partition_sequence.begin(), partition_sequence.end(), // Values to search 
    cut_row_id_to_cell_offset_.begin() // Index?
  );

  // copy last value (total num_cells) to offset vector
  cudaMemcpy(
    cut_row_id_to_cell_offset_.data() + cut_row_id_to_cell_offset_.size() - 1, 
    &num_cells, 
    sizeof(int), 
    cudaMemcpyHostToDevice);

  cut_row_id_to_num_cells_.resize(num_cut_rows_);
  thrust::adjacent_difference( // thrust::adjacent_difference -> output[i] = x[i] - x[i-1]
    cut_row_id_to_cell_offset_.begin() + 1, cut_row_id_to_cell_offset_.begin() + num_cut_rows_ + 1,
    cut_row_id_to_num_cells_.begin());

  //printf("NumCutRow: %d\n", num_cut_rows_);
  //printVector(cut_row_id_to_cell_offset_, "Cell Offset");
  //printVector(cut_row_id_to_num_cells_, "NumCell");
}

void
CudaPostProcess::sortCellsEachCutRow()
{
  const auto& cell_w_in_grid = database_->getCellWidthInGrid();
  const auto& cell_h_in_grid = database_->getCellHeightInGrid();
  const auto& cell_lx_in_grid = database_->getCellLxInGrid();
  const auto& cell_ly_in_grid = database_->getCellLyInGrid();

  int max_cell_in_cut_row = computeVectorMax(cut_row_id_to_num_cells_);
  int max_cell_in_cut_row_padded = findMinPowerOfTwoNoLessThan(max_cell_in_cut_row);

  if(max_cell_in_cut_row_padded > 12288)
    printf("MaxCell: %d MaxCellPadded: %d\n", max_cell_in_cut_row, max_cell_in_cut_row_padded);
  assert(max_cell_in_cut_row_padded <= 12288);

  const size_t shared_mem_size = 2 * sizeof(int) * max_cell_in_cut_row_padded;
  const int num_threads_per_block = 16;

  sortCellsEachCutRowKernel<<<num_cut_rows_, num_threads_per_block, shared_mem_size>>>(
    max_cell_in_cut_row_padded, /* buffer_size */
    cut_row_id_to_num_cells_.data(),
    cut_row_id_to_cell_offset_.data(),
    cell_lx_in_grid.data(),
    cell_grouped_by_cut_row_.data());

  /*
  debugSortKernel<<<1, 1>>>(
    num_cut_rows_,
    cut_row_id_to_num_cells_.data(),
    cut_row_id_to_cell_offset_.data(),
    cut_row_id_to_grid_y_.data(),
    cut_row_id_to_grid_x_bgn_.data(),
    cut_row_id_to_grid_x_end_.data(),
    cell_lx_in_grid.data(),
    cell_ly_in_grid.data(),
    cell_w_in_grid.data(),
    cell_grouped_by_cut_row_.data());
    */
}

void
CudaPostProcess::solveQP()
{
  qp_solver::CudaQpSolver solver(params_, database_);

  solver.initializeSolver(
    cut_row_id_to_grid_x_bgn_,
    cut_row_id_to_grid_x_end_,
    cut_row_id_to_cell_offset_,
    cut_row_id_to_num_cells_,
    cell_grouped_by_cut_row_);

  solver.solve();
}

}
