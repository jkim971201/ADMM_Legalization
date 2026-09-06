#include "legalizer/Legalizer.h"
#include "objects/LGCell.h"
#include "objects/LGCellTechInfo.h"
#include "Grid.h"
#include "GlobalUtil.h"
#include "LGHyperParameters.h"

#include "util/Chrono.h"

#include <cassert>
#include <cmath>
#include <queue>
#include <unordered_set>

namespace legalizer
{

void
Legalizer::initialLegalize()
{
  auto preplace_start = util::getChronoNow();

  const int num_cells = static_cast<int>(cells_.size());
  #pragma omp parallel for num_threads(params_->num_cpu_threads)
  for(int i = 0; i < num_cells; i++)
    preplaceCell(cells_[i]);

  // Run projection for hopeless cells
  #pragma omp parallel for num_threads(params_->num_cpu_threads)
  for(int i = 0; i < num_cells; i++)
  {
    LGCell& cell = cells_[i];
    if(checkHopeless(cell) == true)
    {
      int cell_group_index = cell.getGroupIndex();
      auto& regions_same_group = index2regions_[cell_group_index];
    
      std::vector<ProjectionResult> candidate_pts;
      for(const auto region : regions_same_group)
        candidate_pts.push_back(findProjectionPoint(cell, region.get()));

      const auto [best_lx, best_ly, best_disp] = *std::min_element(candidate_pts.begin(), candidate_pts.end(),
        [] (const ProjectionResult& lhs, const ProjectionResult& rhs) { return lhs.disp < rhs.disp; });

      if(best_disp < k_int_max)
      {
        cell.setLx(best_lx);
        cell.setLy(best_ly);
      }
      else
      {
        printf("No available space for cell %s\n", cell.getName().data());
        assert(0);
      }
    }
  }

  std::string_view max_name;
  int max_disp_dbu = -1;
  for(const auto& cell : cells_)
  {
    int disp = cell.getDisp();
    if(max_disp_dbu < disp)
    {
      max_disp_dbu = disp;
      max_name = cell.getName();
    }
  }
  printf("MaxDisp: %5.2f (sites) Cell: %s\n", double(max_disp_dbu) / double(site_width_), max_name.data());

  printf("preplace           finished (takes %5.2f s)\n", util::evalTime(preplace_start));
}

void
Legalizer::preplaceCell(LGCell& cell)
{
  roundToNearsetPoint(cell);
  const auto& [status, region] = grid_->checkPixelStatus(cell);
  if(status == PixelStatus::VALID)
    return;
  else if(status == PixelStatus::BLOCKAGE or status == PixelStatus::FENCE)
    moveOutsideUncorrectRegion(cell, region.get());
}

void
Legalizer::moveOutsideUncorrectRegion(LGCell& cell, const LGRect* rect)
{
  std::unordered_set<int> visited_grid_pts;
  const int y_size = grid_->getYSize();
  const int cell_group_index = cell.getGroupIndex();

  std::vector<ProjectionResult> candidate_pts;
  const auto& regions = index2regions_[cell_group_index];
  for(const auto region : regions)
    candidate_pts.push_back(findProjectionPoint(cell, region.get()));

  const auto [best_lx, best_ly, best_disp] = *std::min_element(candidate_pts.begin(), candidate_pts.end(),
    [] (const ProjectionResult& lhs, const ProjectionResult& rhs) { return lhs.disp < rhs.disp; });
    
  if(best_disp < k_int_max)
  {
    cell.setLx(best_lx);
    cell.setLy(best_ly);
  }

  const auto& [_stat, _region] = grid_->checkPixelStatus(cell);
  if(_stat != PixelStatus::VALID)
    printf("Preplace fail: %s group_index: %d\n", cell.getName().data(), cell_group_index);
}

void
Legalizer::roundToNearsetPoint(LGCell& cell)
{
  const int lx = cell.getLx();
  const int ly = cell.getLy();

  const int cell_w = cell.getWidth();
  const int cell_w_in_grid = cell_w / site_width_;

  const int cell_h = cell.getHeight();
  const int cell_h_in_grid = cell_h / row_height_;

  const int x_size = grid_->getXSize();
  const int y_size = grid_->getYSize();

  // For X
  int grid_x_min = std::min(x_size - cell_w_in_grid, std::max(0, grid_->getGridStartX(lx)));
  int dbu_x1 = grid_->getDbuX(grid_x_min);
  int dbu_x2 = grid_->getDbuX(std::min(x_size - cell_w_in_grid, grid_x_min + 1));
  int new_lx = (std::abs(lx - dbu_x1) < std::abs(lx - dbu_x2)) ? dbu_x1 : dbu_x2;

  // For Y
  int grid_y_min = std::min(y_size - cell_h_in_grid, std::max(0, grid_->getGridStartY(ly)));
  int dbu_y1 = grid_->getDbuY(grid_y_min);
  int dbu_y2 = grid_->getDbuY(std::min(y_size - cell_h_in_grid, grid_y_min + 1));
  if(dbu_y1 == dbu_y2) // This can occur when grid_y_min is the top (or N-highest when N x Height cell) row.
    dbu_y2 = grid_->getDbuY(std::max(0, grid_y_min - 1));
  
  int new_ly = (std::abs(ly - dbu_y1) < std::abs(ly - dbu_y2)) ? dbu_y1 : dbu_y2;
  if(cell_h_in_grid % 2 == 0)
  {
    const auto& grid_y_to_is_vdd_up = grid_->getGridYToIsVddUp();
    const auto cell_tech_info = cell.getTechInfo();
    const int is_vss_bottom = cell_tech_info->hasGroundAtBottom() == true ? 1 : 0;
    const int new_ly_in_grid = grid_->getGridStartY(new_ly);
    if(grid_y_to_is_vdd_up[new_ly_in_grid] != is_vss_bottom)
    {
      if(new_ly == dbu_y1)
        new_ly = (dbu_y2 + cell_h <= core_uy_) ? dbu_y2 : dbu_y1 - row_height_;
      else if(new_ly == dbu_y2)
        new_ly = (dbu_y1 >= core_ly_) ? dbu_y1 : dbu_y2 + row_height_;
    }
  }

  cell.setLx(new_lx);
  cell.setLy(new_ly);
}

bool 
Legalizer::checkHopeless(LGCell& cell) 
{
  const int num_directions = static_cast<int>(direction_default_.size());
  const auto& grid_y_to_is_vdd_up = grid_->getGridYToIsVddUp();

  const auto cell_tech_info = cell.getTechInfo();
  const int is_vss_bottom = cell_tech_info->hasGroundAtBottom() == true ? 1 : 0;

  const int cell_h_dbu = cell.getHeight();
  const int cell_h_grid = cell_h_dbu / row_height_;

  const int cur_lx_dbu = cell.getLx();
  const int cur_ly_dbu = cell.getLy();

  const int cur_lx_grid = grid_->getGridStartX(cur_lx_dbu);
  const int cur_ly_grid = grid_->getGridStartY(cur_ly_dbu);

  bool is_hopeless = true;
  for(int j = 0; j < num_directions; j++)
  {
    int delta_x = direction_default_[2 * j + 0];
    int delta_y = direction_default_[2 * j + 1];

    int new_lx_grid = cur_lx_grid + delta_x;
    int new_ly_grid = cur_ly_grid + delta_y;

    int new_lx_dbu = grid_->getDbuX(new_lx_grid);
    int new_ly_dbu = grid_->getDbuY(new_ly_grid);

    cell.setLx(new_lx_dbu);
    cell.setLy(new_ly_dbu);

    if(isOutOfCore(cell) == true)
      continue;

    const auto [status, region] = grid_->checkPixelStatus(cell);
    if(grid_y_to_is_vdd_up[new_ly_grid] == is_vss_bottom and status == PixelStatus::VALID)
    {
      is_hopeless = false;
      break;
    }
  }

  cell.setLx(cur_lx_dbu);
  cell.setLy(cur_ly_dbu);

  return is_hopeless;
}

bool 
Legalizer::checkOverlap(const LGRect& rect1, const LGRect& rect2) const
{
  const int lx1 = rect1.getLx();
  const int ly1 = rect1.getLy();
  const int ux1 = rect1.getUx();
  const int uy1 = rect1.getUy();

  const int lx2 = rect2.getLx();
  const int ly2 = rect2.getLy();
  const int ux2 = rect2.getUx();
  const int uy2 = rect2.getUy();

  if(lx2 >= ux1) return false; 
  if(lx1 >= ux2) return false; 
  if(ly2 >= uy1) return false; 
  if(ly1 >= uy2) return false; 
  return true;
}

bool
Legalizer::isOutOfCore(const LGRect& rect) const
{
  const int lx = rect.getLx();
  const int ly = rect.getLy();
  const int ux = rect.getUx();
  const int uy = rect.getUy();

  if(ux > core_ux_) return true;
  if(uy > core_uy_) return true;
  if(lx < core_lx_) return true;
  if(ly < core_ly_) return true;
  return false;
}

bool
Legalizer::isOutOfRegion(const LGRect* region, const LGCell& cell, int cell_lx_dbu, int cell_ly_dbu) const
{
  const int region_lx_dbu = region->getLx();
  const int region_ly_dbu = region->getLy();
  const int region_ux_dbu = region->getUx();
  const int region_uy_dbu = region->getUy();

  const int cell_w_dbu = cell.getWidth();
  const int cell_h_dbu = cell.getHeight();

  bool lx_check = cell_lx_dbu >= region_lx_dbu ? true : false;
  bool ly_check = cell_ly_dbu >= region_ly_dbu ? true : false;

  bool ux_check = (cell_lx_dbu + cell_w_dbu) <= region_ux_dbu ? true : false;
  bool uy_check = (cell_ly_dbu + cell_h_dbu) <= region_uy_dbu ? true : false;

  if(lx_check and ly_check and ux_check and uy_check)
    return false;
  else
    return true;
}

ProjectionResult
Legalizer::findProjectionPoint(const LGCell& cell, const LGRect* region) const
{
  const int cell_lx = cell.getLx();
  const int cell_ly = cell.getLy();
  const int cell_ux = cell.getUx();
  const int cell_uy = cell.getUy();
  const int cell_w  = cell.getWidth();
  const int cell_h  = cell.getHeight();

  const int cell_w_in_grid = cell_w / site_width_;
  const int cell_h_in_grid = cell_h / row_height_;

  const int cell_init_lx = cell.getInitLx();
  const int cell_init_ly = cell.getInitLy();

  const auto [grid_x1, grid_y1, grid_x2, grid_y2] = grid_->getEndPointIndexShrunk(region);
  if(grid_y1 > grid_y2 or grid_x1 > grid_x2)
    return {-1, -1, k_int_max};

  const int region_lx = grid_->getDbuX(grid_x1);
  const int region_ly = grid_->getDbuY(grid_y1);
  const int region_ux = grid_->getDbuX(grid_x2);
  const int region_uy = grid_->getDbuY(grid_y2);

  int new_lx = cell_lx;
  int new_ly = cell_ly;

  if(cell_ux > region_ux)
    new_lx = grid_->getDbuX(grid_x2 - cell_w_in_grid + 1);
  else if(cell_lx < region_lx)
    new_lx = grid_->getDbuX(grid_x1);

  if(cell_uy > region_uy)
    new_ly = grid_->getDbuY(grid_y2 - cell_h_in_grid + 1);
  else if(cell_ly < region_ly)
    new_ly = grid_->getDbuY(grid_y1);

  if(cell_h_in_grid % 2 == 0)
  {
    const int new_ly_in_grid = grid_->getGridStartY(new_ly);
    const auto cell_tech_info = cell.getTechInfo();
    const int is_vss_bottom = cell_tech_info->hasGroundAtBottom() == true ? 1 : 0;
    const auto& grid_y_to_is_vdd_up = grid_->getGridYToIsVddUp();
    if(grid_y_to_is_vdd_up[new_ly_in_grid] != is_vss_bottom)
    {
      std::array<int, 2> dbu_y_candidates = {new_ly - row_height_, new_ly + row_height_};
      int best_disp = std::numeric_limits<int>::max();
      int best_ly = -1;

      for(int dbu_y_cand : dbu_y_candidates)
      {
        int disp_temp = std::abs(new_lx - cell_init_lx) + std::abs(dbu_y_cand - cell_init_ly);
        if(isOutOfRegion(region, cell, new_lx, dbu_y_cand) == false and disp_temp < best_disp)
        {
          best_disp = disp_temp;
          best_ly = dbu_y_cand;
        }
      }

      new_ly = best_ly;
    }
  }

  bool is_invalid = false;
  bool is_out_of_region = isOutOfRegion(region, cell, new_lx, new_ly);

  // this is for superblue11
  if(cell_h_in_grid % 2 == 0 and is_out_of_region == false)
  {
    const auto& pixels = grid_->getPixels();
    const int new_x_in_grid = grid_->getGridStartX(new_lx);
    const int new_y_in_grid = grid_->getGridStartY(new_ly);
    //printf("x_min: %d y_min: %d x_max: %d y_max: %d\n", 
    //  new_x_in_grid, new_y_in_grid, new_x_in_grid + cell_w_in_grid - 1, new_y_in_grid + cell_h_in_grid - 1);
    for(int x = new_x_in_grid; x < new_x_in_grid + cell_w_in_grid; x++)
      for(int y = new_y_in_grid; y < new_y_in_grid + cell_h_in_grid; y++)
        if(pixels[x][y].valid == false)
          is_invalid = true;
  }

  int final_disp = std::abs(cell_init_lx - new_lx) + std::abs(cell_init_ly - new_ly);
  if(is_out_of_region == false and is_invalid == false)
    return {new_lx, new_ly, final_disp};
  else
    return {-1, -1, k_int_max};
}

bool
Legalizer::checkRowViolation(const LGCell& cell) const
{
  const int cell_ly_dbu = cell.getLy();
  const int grid_y = grid_->getGridStartY(cell_ly_dbu);
  const int cell_h_in_grid = cell.getHeight() / row_height_;
  const auto cell_tech_info = cell.getTechInfo();
  const int is_vss_bottom = cell_tech_info->hasGroundAtBottom() == true ? 1 : 0;
  const auto& grid_y_to_is_vdd_up = grid_->getGridYToIsVddUp();
  if(cell_h_in_grid % 2 == 0 and grid_y_to_is_vdd_up[grid_y] != is_vss_bottom)
    return true;
  else
    return false;
}

}
