#include "CudaDatabase.h"
#include "DatabaseKernels.hpp"

#include "Grid.h"
#include "LGHyperParameters.h"
#include "objects/LGVioStat.h"
#include "objects/LGCell.h"
#include "objects/LGCellTechInfo.h"

#include <fstream>

#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>

#include "cuda_linalg/CudaVectorAlgebra.h"

namespace legalizer
{

CudaDatabase::CudaDatabase(
  const std::vector<std::shared_ptr<LGRegion>>& regions,
  const std::shared_ptr<LGHyperParameters> param, 
  const std::shared_ptr<Grid> grid) : param_(param), grid_(grid)
{
  x_size_ = grid_->getXSize();
  y_size_ = grid_->getYSize();
  
  site_width_ = grid_->getSiteWidth();
  row_height_ = grid_->getRowHeight();

  core_lx_ = grid_->getLx();
  core_ly_ = grid_->getLy();

  copyGridInfo(regions);

  diagnosePixels();
}

void
CudaDatabase::copyFromHostCells(
  const std::vector<LGCell>& cells,
  const std::vector<std::shared_ptr<LGCellTechInfo>>& tech_infos)
{
  num_cells_ = cells.size();

  std::vector<int> host_cell_vector(cells.size());

  // Single-thread is faster than multi-thread
  // (because of memory access pattern?)

  /* Cell Original Lx and Ly (dbu) */
  for(int i = 0; i < num_cells_; i++) 
    host_cell_vector[i] = cells[i].getInitLx();
  cell_id_to_original_lx_in_dbu_ = host_cell_vector;

  for(int i = 0; i < num_cells_; i++) 
    host_cell_vector[i] = cells[i].getInitLy();
  cell_id_to_original_ly_in_dbu_ = host_cell_vector;

  /* Cell Lx and Ly (dbu) */
  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getLx();
  cell_id_to_lx_in_dbu_ = host_cell_vector;

  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getLy();
  cell_id_to_ly_in_dbu_ = host_cell_vector;

  /* Cell W and H (grid) */
  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getWidth() / site_width_;
  cell_id_to_width_in_grid_ = host_cell_vector;

  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getHeight() / row_height_;
  cell_id_to_height_in_grid_ = host_cell_vector;

  /* Cell Group Index */
  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getGroupIndex();
  cell_id_to_group_index_ = host_cell_vector;

  /* Cell Orient */
  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getTechInfo()->hasGroundAtBottom() == true ? 1 : 0;
  cell_id_to_has_ground_at_bottom_ = host_cell_vector;

  /* Cell Edge Type */
  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getLEdgeType();
  cell_id_to_ledge_type_ = host_cell_vector;

  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getREdgeType();
  cell_id_to_redge_type_ = host_cell_vector;

  /* Cell to Number of Pin Short/Access, Edge Spacing Violations */
  cell_id_to_num_short_vio_.resize(cells.size());
  cell_id_to_num_access_vio_.resize(cells.size());
  cell_id_to_num_edge_spacing_vio_.resize(cells.size());

  /* Cell Lx and Ly (grid) */
  cell_id_to_lx_in_grid_.resize(cells.size());
  cell_id_to_ly_in_grid_.resize(cells.size());

  convertDbuToGrid(
    cell_id_to_lx_in_dbu_,
    cell_id_to_ly_in_dbu_,
    cell_id_to_lx_in_grid_,
    cell_id_to_ly_in_grid_);

  /* Cell Orientation-related */
  cell_id_to_need_rotation_.resize(cells.size());
  cell_id_to_need_flip_.resize(cells.size());

  // Technology Info
  std::vector<int> host_pin_offset;
  std::vector<int> host_pin_bgn;
  std::vector<int> host_pin_end;
  std::vector<int> host_pin_layer;

  int cell_pin_id = 0;
  for(const auto& cell_tech_info : tech_infos)
  {
    const auto& pin_shapes = cell_tech_info->getPinShapes();
    host_pin_offset.push_back(cell_pin_id);
    for(const auto& [pin_layer, pin_bgn, pin_end] : pin_shapes)
    {
      host_pin_bgn.push_back(pin_bgn);
      host_pin_end.push_back(pin_end);
      host_pin_layer.push_back(pin_layer);
      cell_pin_id++;
    }
  }
  host_pin_offset.push_back(cell_pin_id);

  for(int i = 0; i < num_cells_; i++)
    host_cell_vector[i] = cells[i].getTechInfo()->getID();

  cell_id_to_cell_macro_id_ = host_cell_vector;
  cell_macro_id_to_pin_offset_ = host_pin_offset;
  pin_id_to_bgn_   = host_pin_bgn;
  pin_id_to_end_   = host_pin_end;
  pin_id_to_layer_ = host_pin_layer;
}

void
CudaDatabase::copyGridInfo(const std::vector<std::shared_ptr<LGRegion>>& regions)
{
  const int num_pixels = x_size_ * y_size_;
  pixel_id_to_usage_.resize(num_pixels);

  std::vector<int> host_pixel_vector(num_pixels);
  
  const auto& pixels = grid_->getPixels();

  /* pixel_id to is_valid */
  for(int x = 0; x < x_size_; x++)
    for(int y = 0; y < y_size_; y++)
      host_pixel_vector[x * y_size_ + y] = pixels[x][y].valid == true ? 1 : 0;
  pixel_id_to_valid_ = host_pixel_vector;

  /* pixel_id to group_id */
  for(int x = 0; x < x_size_; x++)
    for(int y = 0; y < y_size_; y++)
      host_pixel_vector[x * y_size_ + y] = pixels[x][y].group_index;
  pixel_id_to_group_id_ = host_pixel_vector;

  /* pixel_id to power metal layer */
  for(int x = 0; x < x_size_; x++)
    for(int y = 0; y < y_size_; y++)
      host_pixel_vector[x * y_size_ + y] = pixels[x][y].power_metal_layer;
  pixel_id_to_power_metal_layer_ = host_pixel_vector;

  /* pixel_id to is isolated region */
  for(int x = 0; x < x_size_; x++)
  {
    for(int y = 0; y < y_size_; y++)
    {
      auto region = pixels[x][y].region;
      if(region != nullptr)
        host_pixel_vector[x * y_size_ + y] = region->getIndex();
      else
        host_pixel_vector[x * y_size_ + y] = -1;
    }
  }

  pixel_id_to_region_id_ = host_pixel_vector;

  /* region information */
  const int num_regions = static_cast<int>(regions.size());

  std::vector<int> host_region_grid_x_min(num_regions);
  std::vector<int> host_region_grid_y_min(num_regions);
  std::vector<int> host_region_grid_x_max(num_regions);
  std::vector<int> host_region_grid_y_max(num_regions);

  for(int i = 0; i < num_regions; i++)
  {
    const auto region = regions[i];
    const auto [x_min, y_min, x_max, y_max] = grid_->getEndPointIndexShrunk(region.get());
    host_region_grid_x_min[i] = x_min;
    host_region_grid_y_min[i] = y_min;
    host_region_grid_x_max[i] = x_max;
    host_region_grid_y_max[i] = y_max;
  }

  region_id_to_grid_x_min_ = host_region_grid_x_min;
  region_id_to_grid_y_min_ = host_region_grid_y_min;
  region_id_to_grid_x_max_ = host_region_grid_x_max;
  region_id_to_grid_y_max_ = host_region_grid_y_max;

  /* pixel_id to num edge type */
  pixel_id_to_num_ledge_type1_.resize(num_pixels);
  pixel_id_to_num_ledge_type2_.resize(num_pixels);
  pixel_id_to_num_redge_type1_.resize(num_pixels);
  pixel_id_to_num_redge_type2_.resize(num_pixels);

  /* pixel_id_to_shape */
  pixel_id_to_shape_.resize(num_pixels);
  // will be initialized by diagnosePixels()

  /* grid_y to is_vdd_up */
  grid_y_to_is_vdd_up_ = grid_->getGridYToIsVddUp();
  //printVector(grid_y_to_is_vdd_up_, "vdd_up");
}

void
CudaDatabase::copyEdgeSpacingRules(const std::vector<int>& edge_spacing_rules)
{
  edge_spacing_rules_ = edge_spacing_rules;
}

void
CudaDatabase::convertDbuToGrid(
  const cuda_linalg::CudaVector<int>& lx_in_dbu,
  const cuda_linalg::CudaVector<int>& ly_in_dbu,
        cuda_linalg::CudaVector<int>& lx_in_grid,
        cuda_linalg::CudaVector<int>& ly_in_grid)
{
  thrust::transform(
    lx_in_dbu.begin(), lx_in_dbu.end(), 
    lx_in_grid.begin(), 
    dbu_to_grid_functor(core_lx_, site_width_));

  thrust::transform(
    ly_in_dbu.begin(), ly_in_dbu.end(), 
    ly_in_grid.begin(), 
    dbu_to_grid_functor(core_ly_, row_height_));
}

void
CudaDatabase::convertGridToDbu(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
        cuda_linalg::CudaVector<int>& lx_in_dbu,
        cuda_linalg::CudaVector<int>& ly_in_dbu)
{
  thrust::transform(
    lx_in_grid.begin(), lx_in_grid.end(), 
    lx_in_dbu.begin(), 
    grid_to_dbu_functor(core_lx_, site_width_));

  thrust::transform(
    ly_in_grid.begin(), ly_in_grid.end(), 
    ly_in_dbu.begin(), 
    grid_to_dbu_functor(core_ly_, row_height_));
}

void
CudaDatabase::copyToHostCells(std::vector<LGCell>& cells)
{
  std::vector<int> host_rotation(num_cells_);

  std::vector<int> host_cell_x_in_grid(num_cells_);
  std::vector<int> host_cell_y_in_grid(num_cells_);

  thrust::copy(
    cell_id_to_lx_in_grid_.begin(), cell_id_to_lx_in_grid_.end(),
    host_cell_x_in_grid.begin());

  thrust::copy(
    cell_id_to_ly_in_grid_.begin(), cell_id_to_ly_in_grid_.end(),
    host_cell_y_in_grid.begin());

  thrust::copy(
    cell_id_to_need_rotation_.begin(), cell_id_to_need_rotation_.end(),
    host_rotation.begin());

  for(int i = 0; i < num_cells_; i++)
  {
    int grid_x = host_cell_x_in_grid[i];
    int grid_y = host_cell_y_in_grid[i];

    int dbu_x = grid_->getDbuX(grid_x);
    int dbu_y = grid_->getDbuY(grid_y);

    int need_rotation = host_rotation[i];

    cells[i].setLx(dbu_x);
    cells[i].setLy(dbu_y);
    if(need_rotation == 1)
      cells[i].setOrient(Orient::FS);
  }
}

void 
CudaDatabase::commitPlaceInGrid(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid)
{
  thrust::copy(lx_in_grid.begin(), lx_in_grid.end(), cell_id_to_lx_in_grid_.begin());
  thrust::copy(ly_in_grid.begin(), ly_in_grid.end(), cell_id_to_ly_in_grid_.begin());

  convertGridToDbu(
    cell_id_to_lx_in_grid_,
    cell_id_to_ly_in_grid_,
    cell_id_to_lx_in_dbu_,
    cell_id_to_ly_in_dbu_);
}

LGVioStat
CudaDatabase::checkPlace()
{
  //printf("=== CheckPlace ===\n");
  // NOTE : We assume both lx_in_dbu, ly_in_dbu and lx_in_grid, ly_in_grid
  // are synchronized with each other.
  
  /* 1. Check Overlap */
  updateUsageInGrid(cell_id_to_lx_in_grid_, cell_id_to_ly_in_grid_);
  int ovf_sum = computeOverflow();

  /* 2. Row Orientation */
  cuda_linalg::CudaVector<int> cell_id_to_num_row_orient_vio(num_cells_);
  countRowOrientVio(cell_id_to_lx_in_grid_, cell_id_to_ly_in_grid_, cell_id_to_num_edge_spacing_vio_);
  int num_row_orient_vio = computeVectorSum(cell_id_to_num_edge_spacing_vio_);

  /* 3. Pin Short */
  countPinVio(cell_id_to_lx_in_grid_, cell_id_to_ly_in_grid_, cell_id_to_num_short_vio_, cell_id_to_num_access_vio_);
  int num_pin_short_vio = computeVectorSum(cell_id_to_num_short_vio_);
  int num_pin_access_vio = computeVectorSum(cell_id_to_num_access_vio_);

  /* 4. Edge Spacing */
  countEdgeSpacingVio(cell_id_to_lx_in_grid_, cell_id_to_ly_in_grid_, cell_id_to_num_edge_spacing_vio_);
  int num_edge_spacing_vio = computeVectorSum(cell_id_to_num_edge_spacing_vio_) / 2;
  // Like NBLG, we should divide by 2 because we don't want to double count the same pairs.

  return {ovf_sum, num_row_orient_vio, num_pin_short_vio, num_pin_access_vio, num_edge_spacing_vio};
}

void
CudaDatabase::countPinVio(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
        cuda_linalg::CudaVector<int>& cell_id_to_num_pin_short_vio,
        cuda_linalg::CudaVector<int>& cell_id_to_num_pin_access_vio)
{
  cell_id_to_num_pin_short_vio.fillZero();
  cell_id_to_num_pin_access_vio.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  countPinVioKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    x_size_,
    y_size_,
    lx_in_grid.data(),
    ly_in_grid.data(),
    cell_id_to_width_in_grid_.data(),
    cell_id_to_height_in_grid_.data(),
    cell_id_to_cell_macro_id_.data(),
    cell_macro_id_to_pin_offset_.data(),
    cell_id_to_has_ground_at_bottom_.data(),
    pin_id_to_bgn_.data(),
    pin_id_to_end_.data(),
    pin_id_to_layer_.data(),
    pixel_id_to_power_metal_layer_.data(),
    grid_y_to_is_vdd_up_.data(),
    cell_id_to_num_pin_short_vio.data(),
    cell_id_to_num_pin_access_vio.data());
}

void
CudaDatabase::countEdgeSpacingVio(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
        cuda_linalg::CudaVector<int>& cell_id_to_num_edge_spacing_vio)
{
  cell_id_to_num_edge_spacing_vio.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  countEdgeSpacingVioKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    x_size_,
    y_size_,
    lx_in_grid.data(),
    ly_in_grid.data(),
    cell_id_to_width_in_grid_.data(),
    cell_id_to_height_in_grid_.data(),
    cell_id_to_ledge_type_.data(),
    cell_id_to_redge_type_.data(),
    cell_id_to_has_ground_at_bottom_.data(),
    grid_y_to_is_vdd_up_.data(),
    pixel_id_to_num_ledge_type1_.data(),
    pixel_id_to_num_ledge_type2_.data(),
    pixel_id_to_num_redge_type1_.data(),
    pixel_id_to_num_redge_type2_.data(),
    edge_spacing_rules_.data(),
    cell_id_to_num_edge_spacing_vio.data());

  //printf("Sum1: %d\n", computeVectorSum(pixel_id_to_num_ledge_type1_));
  //printf("Sum2: %d\n", computeVectorSum(pixel_id_to_num_redge_type1_));
}

void
CudaDatabase::countRowOrientVio(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid,
        cuda_linalg::CudaVector<int>& cell_id_to_num_row_orient_vio)
{
  cell_id_to_num_row_orient_vio.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  countRowOrientVioKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    ly_in_grid.data(),
    cell_id_to_height_in_grid_.data(),
    cell_id_to_has_ground_at_bottom_.data(),
    cell_id_to_need_rotation_.data(),
    grid_y_to_is_vdd_up_.data(),
    cell_id_to_num_row_orient_vio.data());
}

void
CudaDatabase::updateUsageInGrid(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid)
{
  pixel_id_to_usage_.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  computeUsageGridKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    x_size_,
    y_size_,
    cell_id_to_width_in_grid_.data(),
    cell_id_to_height_in_grid_.data(),
    lx_in_grid.data(),
    ly_in_grid.data(),
    pixel_id_to_usage_.data());
}

void
CudaDatabase::updatePixelEdgeMapInGrid(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid)
{
  pixel_id_to_num_ledge_type1_.fillZero();
  pixel_id_to_num_ledge_type2_.fillZero();
  pixel_id_to_num_redge_type1_.fillZero();
  pixel_id_to_num_redge_type2_.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  computePixelEdgeMapKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    x_size_,
    y_size_,
    cell_id_to_width_in_grid_.data(),
    cell_id_to_height_in_grid_.data(),
    lx_in_grid.data(),
    ly_in_grid.data(),
    cell_id_to_ledge_type_.data(),
    cell_id_to_redge_type_.data(),
    cell_id_to_has_ground_at_bottom_.data(),
    grid_y_to_is_vdd_up_.data(),
    pixel_id_to_num_ledge_type1_.data(),
    pixel_id_to_num_ledge_type2_.data(),
    pixel_id_to_num_redge_type1_.data(),
    pixel_id_to_num_redge_type2_.data());

  /*
  std::vector<int> host_vec(pixel_id_to_num_ledge_type1_.size(), 0);

  thrust::copy(
    pixel_id_to_num_redge_type1_.begin(), 
    pixel_id_to_num_redge_type1_.end(),
    host_vec.begin());

  int sum1 = 0;
  std::ofstream output1("edge_log1.txt");
  for(int x = 0; x < x_size_; x++)
  {
    for(int y = 0; y < y_size_; y++)
    {
      int val = host_vec[x * y_size_ + y];
      if(val >= 1)
        output1 << x << " " << y << std::endl;
      sum1 += val;
    }
  }

  thrust::copy(
    pixel_id_to_num_ledge_type1_.begin(), 
    pixel_id_to_num_ledge_type1_.end(),
    host_vec.begin());

  int sum2 = 0;
  std::ofstream output2("edge_log2.txt");
  for(int x = 0; x < x_size_; x++)
  {
    for(int y = 0; y < y_size_; y++)
    {
      int val = host_vec[x * y_size_ + y];
      if(val >= 1)
        output2 << x << " " << y << std::endl;
      sum2 += val;
    }
  }

  printf("Sum1: %d Sum2: %d\n", sum1, sum2);
  */
}

void
CudaDatabase::updateUsageInDbu(
  const cuda_linalg::CudaVector<int>& lx_in_dbu,
  const cuda_linalg::CudaVector<int>& ly_in_dbu)
{
  pixel_id_to_usage_.fillZero();

  int num_thread_per_block = 512;
  int num_block = (num_cells_ + num_thread_per_block - 1) / num_thread_per_block;

  computeUsageDbuKernel<<<num_block, num_thread_per_block>>>(
    num_cells_,
    core_lx_,
    core_ly_,
    site_width_,
    row_height_,
    x_size_,
    y_size_,
    cell_id_to_width_in_grid_.data(),
    cell_id_to_height_in_grid_.data(),
    lx_in_dbu.data(),
    ly_in_dbu.data(),
    pixel_id_to_usage_.data());
}

int
CudaDatabase::computeOverflow() const
{
  int ovf_sum = thrust::transform_reduce(
    pixel_id_to_usage_.begin(), pixel_id_to_usage_.end(),
    overflow_functor(),
    0,
    thrust::plus<int>());

  return ovf_sum;
}

float
CudaDatabase::computeDisplace() const
{
  return this->computeDisplace(cell_id_to_lx_in_grid_, cell_id_to_ly_in_grid_);
}

float
CudaDatabase::computeDisplace(
  const cuda_linalg::CudaVector<int>& lx_in_grid,
  const cuda_linalg::CudaVector<int>& ly_in_grid) const
{
  auto zip_bgn = thrust::make_zip_iterator(
    thrust::make_tuple(
      lx_in_grid.begin(),
      ly_in_grid.begin(),
      cell_id_to_original_lx_in_dbu_.begin(),
      cell_id_to_original_ly_in_dbu_.begin()));

  auto zip_end = thrust::make_zip_iterator(
    thrust::make_tuple(
      lx_in_grid.end(),
      ly_in_grid.end(),
      cell_id_to_original_lx_in_dbu_.end(),
      cell_id_to_original_ly_in_dbu_.end()));

  float disp = thrust::transform_reduce(zip_bgn, zip_end,
    displace_functor(core_lx_, core_ly_, site_width_, row_height_), 0.0f, thrust::plus<float>());

  return disp;
}

void
CudaDatabase::writeViolationMap(const cuda_linalg::CudaVector<int>& access_vio) const
{
  std::ofstream output("pin_layer.txt");

  std::vector<int> host_vio(access_vio.size());
  thrust::copy(access_vio.begin(), access_vio.end(), host_vio.begin());

  for(int pixel_id = 0; pixel_id < host_vio.size(); pixel_id++)
  {
    int x = pixel_id / y_size_;
    int y = pixel_id % y_size_;
    int vio = host_vio[pixel_id];
    output << x << " " << y << " " << vio << std::endl;
  }
}

void
CudaDatabase::diagnosePixels()
{
  const int num_pixels = x_size_ * y_size_;

  const int num_thread_per_block = 512;

  const int num_block = (num_pixels + num_thread_per_block - 1) / num_thread_per_block;

  constexpr int max_search = 100;
  constexpr int k_ratio = 3;

  diagnosePixelKernel<<<num_block, num_thread_per_block>>>(
    num_pixels,
    x_size_,
    y_size_,
    max_search,
    k_ratio,
    pixel_id_to_valid_.data(),
    pixel_id_to_group_id_.data(),
    pixel_id_to_shape_.data());

  /*
  std::ofstream output("shape_map.txt");

  std::vector<int> host_shape(pixel_id_to_shape_.size());
  thrust::copy(pixel_id_to_shape_.begin(), pixel_id_to_shape_.end(), host_shape.begin());

  for(int pixel_id = 0; pixel_id < host_shape.size(); pixel_id++)
  {
    int x = pixel_id / y_size_;
    int y = pixel_id % y_size_;
    int shape = host_shape[pixel_id];
    output << x << " " << y << " " << shape << std::endl;
  }
  */
}

}
