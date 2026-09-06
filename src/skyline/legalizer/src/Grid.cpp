#include <cmath>
#include <cassert>
#include <numeric>
#include <unordered_set>

#include "db/dbDatabase.h"
#include "db/dbDesign.h"
#include "db/dbNet.h"
#include "db/dbRow.h"
#include "db/dbDie.h"
#include "db/dbRect.h"
#include "db/dbBlockage.h"

#include "Grid.h"

namespace legalizer
{

Grid::Grid(std::shared_ptr<db::dbDatabase> db) : site_width_(0), row_height_(0)
{
  initNonGroupGrid(db);
}

void
Grid::initNonGroupGrid(std::shared_ptr<db::dbDatabase> db)
{
  const auto design = db->getDesign();
  const std::vector<dbInst*>& db_insts = design->getInsts();
  const std::vector<dbRow*>& db_rows   = design->getRows();

  int core_lx = design->coreLx();
  int core_ly = design->coreLy();
  int core_ux = design->coreUx();
  int core_uy = design->coreUy();

  lx_ = core_lx;
  ly_ = core_ly;
  w_ = core_ux - core_lx;
  h_ = core_uy - core_ly;

  std::unordered_set<int> unique_site_width;
  std::unordered_set<int> unique_row_height;
  int min_row_height = std::numeric_limits<int>::max();
  for(const auto& db_row : db_rows)
  {
    min_row_height = std::min(min_row_height, db_row->dy());
    unique_site_width.insert(db_row->siteWidth());
    unique_row_height.insert(db_row->dy());
  }

  assert(unique_site_width.size() == 1);
  assert(unique_row_height.size() == 1);

  row_height_ = min_row_height;
  site_width_ = db_rows.back()->siteWidth(); // TODO: Fix this.

  x_size_ = w_ / site_width_;
  y_size_ = h_ / row_height_;

  grid_y_to_is_vdd_up_.resize(y_size_, 0);
  pixels_.resize(x_size_, std::vector<Pixel>(y_size_, Pixel()));

  for(const auto& db_row : db_rows)
    markDbRow(db_row);

  printf("Grid: (%d,%d)\n", x_size_, y_size_);
}

int
Grid::getGridStartX(int dbu_x) const
{
  return (dbu_x - lx_) / site_width_;
}

int
Grid::getGridStartY(int dbu_y) const
{
  return (dbu_y - ly_) / row_height_;
}

int
Grid::getGridEndX(int dbu_x) const
{
  int quo = (dbu_x - lx_) / site_width_;
  int rem = (dbu_x - lx_) % site_width_;
  int end_x = rem == 0 ? quo - 1 : quo;
  return std::max(std::min(x_size_ - 1, end_x), 0);
}

int
Grid::getGridEndY(int dbu_y) const
{
  int quo = (dbu_y - ly_) / row_height_;
  int rem = (dbu_y - ly_) % row_height_;
  int end_y = rem == 0 ? quo - 1 : quo;
  return std::max(std::min(y_size_ - 1, end_y), 0);
}

std::array<int, 4>
Grid::getEndPointIndexExpanded(const LGRect* rect) const
{
  const int x1 = getGridStartX(rect->getLx());
  const int y1 = getGridStartY(rect->getLy());
  const int x2 = getGridEndX(rect->getUx());
  const int y2 = getGridEndY(rect->getUy());
  return {x1, y1, x2, y2};
}

std::array<int, 4>
Grid::getEndPointIndexShrunk(const LGRect* rect) const
{
  const int x1 = getGridEndX(rect->getLx() + site_width_);
  const int y1 = getGridEndY(rect->getLy() + row_height_);
  const int x2 = getGridStartX(rect->getUx() - site_width_);
  const int y2 = getGridStartY(rect->getUy() - row_height_);

  //if(x2 < x1 or y2 < y1)
  //  printf("rect (%d, %d) - (%d, %d) -> Index (%d, %d) - (%d, %d) CoreLxLy (%d, %d)\n",
  //    rect->getLx(), rect->getLy(), rect->getUx(), rect->getUy(), x1, y1, x2, y2, lx_, ly_);
  //assert(x2 >= x1);
  //assert(y2 >= y1);
  return {x1, y1, x2, y2};
}

int 
Grid::getDbuX(int grid_x) const
{
  return grid_x * site_width_ + lx_;
}

int 
Grid::getDbuY(int grid_y) const
{
  return grid_y * row_height_ + ly_;
}

std::optional<std::pair<int, int>>
Grid::findBestLegalPointNearby(const LGCell& cell, int real_x, int real_y) const
{
  int width_in_site = cell.getWidth() / site_width_;
  int height_in_row = cell.getHeight() / row_height_;

  int initial_x = real_x / site_width_;
  int initial_y = real_y / row_height_;

  std::vector<std::pair<int, int>> grid_pts;

  auto insert_point = [&] (int grid_x, int grid_y)
  {
    if(grid_x >= 0 and grid_x < x_size_ and grid_y >= 0 and grid_y < y_size_)
    {
      if(isLegalMove(cell, getDbuX(grid_x), getDbuY(grid_y)) == true)
        grid_pts.push_back({grid_x, grid_y});
    }
  };

  insert_point(initial_x, initial_y);
  insert_point(initial_x + 1, initial_y);
  insert_point(initial_x - 1, initial_y);
  insert_point(initial_x, initial_y + 1);
  insert_point(initial_x, initial_y - 1);
    
  if(grid_pts.empty() == true)
    return std::nullopt; // return nullopt when fail to find a legal point nearby
  else
  {
    auto compute_disp = [&] (const std::pair<int, int>& move)
    {
      int new_lx = getDbuX(move.first);
      int new_ly = getDbuY(move.second);
      return std::abs(new_lx - cell.getInitLx()) + std::abs(new_ly - cell.getInitLy());
    };

    const auto [best_grid_x, best_grid_y]
      = *std::min_element(grid_pts.begin(), grid_pts.end(), 
        [&] (const auto& pt1, const auto& pt2) { return compute_disp(pt1) < compute_disp(pt2); });

    // convert to real x y coordinates
    std::pair<int, int> best_pt = {getDbuX(best_grid_x), getDbuY(best_grid_y)};
    return best_pt;
  }
};

void 
Grid::markGroupRegion(const std::shared_ptr<LGRegion> region)
{
  const RegionType region_type = region->getType();

  const int region_group_id = region->getGroupIndex();

  const int region_lx = region->getLx();
  const int region_ly = region->getLy();
  const int region_ux = region->getUx();
  const int region_uy = region->getUy();

  const auto [x1_exp, y1_exp, x2_exp, y2_exp] = getEndPointIndexExpanded(region.get());
  const auto [x1, y1, x2, y2] = getEndPointIndexShrunk(region.get());
  // If region width/height is smaller than site_width/row_height, x1/y1 will be larger than x2/y2.

  for(int i = x1_exp; i <= x2_exp; i++)
  {
    for(int j = y1_exp; j <= y2_exp; j++)
    {
      auto& pixel = pixels_[i][j];
      pixel.group_index = region_group_id;
      pixel.h_overlap += site_width_;
      pixel.v_overlap += row_height_;
      pixel.region = region;
      pixel.valid = true;
    }
  }

  //printf("Mark REGION (%d, %d) - (%d, %d)\n", region_lx, region_ly, region_ux, region_uy);

  if(x1_exp == x1 - 1)
  {
    int x_overlap = x1 * site_width_ - region_lx;
    for(int y = y1_exp; y <= y2_exp; y++)
    {
      auto& pixel = pixels_[x1_exp][y];
      pixel.h_overlap += x_overlap - site_width_;
      //printf("  1 Region (%d, %d) - (%d, %d) ", region_lx, region_ly, region_ux, region_uy);
      //printf("(%6d, %3d) h_overlap: %3d -> h_overlap: %d\n", x1_exp, y, x_overlap, pixel.h_overlap);
    }
  }

  if(x2_exp == x2 + 1)
  {
    int x_overlap = region_ux - x2_exp * site_width_;
    for(int y = y1_exp; y <= y2_exp; y++)
    {
      auto& pixel = pixels_[x2_exp][y];
      pixel.h_overlap += x_overlap - site_width_;
      //printf("  2 Region (%d, %d) - (%d, %d) ", region_lx, region_ly, region_ux, region_uy);
      //printf("(%6d, %3d) x_overlap: %3d -> h_overlap: %d\n", x2_exp, y, x_overlap, pixel.h_overlap);
    }
  }

  if(y1_exp == y1 - 1)
  {
    int y_overlap = y1 * row_height_ - region_ly;
    for(int x = x1_exp; x <= x2_exp; x++)
    {
      auto& pixel = pixels_[x][y1_exp];
      pixel.v_overlap += y_overlap - row_height_;
      //printf("  3 Region (%d, %d) - (%d, %d) ", region_lx, region_ly, region_ux, region_uy);
      //printf("(%6d, %3d) y_overlap: %3d -> v_overlap: %d\n", x, y1_exp, y_overlap, pixel.v_overlap);
    }
  }

  if(y2_exp == y2 + 1)
  {
    int y_overlap = region_uy - y2_exp * row_height_;
    for(int x = x1_exp; x <= x2_exp; x++)
    {
      auto& pixel = pixels_[x][y2_exp];
      pixel.v_overlap += y_overlap - row_height_;
      //printf("  4 Region (%d, %d) - (%d, %d) ", region_lx, region_ly, region_ux, region_uy);
      //printf("(%6d, %3d) y_overlap: %3d - > v_overlap: %d\n", x, y2_exp, y_overlap, pixel.v_overlap);
    }
  }
}

void 
Grid::markNonGroupRegion(const std::shared_ptr<LGRegion> region)
{
  const RegionType region_type = region->getType();
  const auto [x1, y1, x2, y2] = getEndPointIndexShrunk(region.get());

  int x1_clamp = std::min(x_size_ - 1, std::max(0, x1));
  int y1_clamp = std::min(y_size_ - 1, std::max(0, y1));
  int x2_clamp = std::min(x_size_ - 1, std::max(0, x2));
  int y2_clamp = std::min(y_size_ - 1, std::max(0, y2));

  // NOTE: Dont change validity here, only link region pointer
  for(int i = x1_clamp; i <= x2_clamp; i++)
  {
    for(int j = y1_clamp; j <= y2_clamp; j++)
    {
      auto& pixel = pixels_[i][j];
      pixel.region = region;
    }
  }
}

void 
Grid::markBlockageRegion(const std::shared_ptr<LGRegion> region)
{
  const RegionType region_type = region->getType();

  const int region_group_id = region->getGroupIndex();

  auto [x1_exp, y1_exp, x2_exp, y2_exp] = getEndPointIndexExpanded(region.get());
  
  if(x1_exp < 0 or x2_exp >= x_size_ or y1_exp < 0 or y2_exp >= y_size_)
  {
    x1_exp = std::max(0, x1_exp);
    y1_exp = std::max(0, y1_exp);
    x2_exp = std::min(x_size_ - 1, x2_exp);
    y2_exp = std::min(y_size_ - 1, y2_exp);
    printf("Warning - region out of core\n");
  }

  for(int i = x1_exp; i <= x2_exp; i++)
  {
    for(int j = y1_exp; j <= y2_exp; j++)
    {
      auto& pixel = pixels_[i][j];
      pixel.valid = false;
      pixel.group_index = region_group_id;
      pixel.region = region;
    }
  }
}

void 
Grid::markDbRow(const db::dbRow* db_row)
{
  LGRect row_rect(db_row->lx(), db_row->ly(), db_row->dx(), db_row->dy());

  const auto [x1, y1, x2, y2] = getEndPointIndexShrunk(&row_rect);

  // If region width/height is smaller than site_width/row_height,
  // x1/y1 will be larger than x2/y2
  if(y1 > y2 or x1 > x2)
    return;

  for(int i = x1; i <= x2; i++)
  {
    for(int j = y1; j <= y2; j++)
    {
      auto& pixel = pixels_[i][j];
      pixel.valid = true;
    }
  }

  const auto row_orient = db_row->orient();
  if(row_orient == Orient::N || row_orient == Orient::FN)
    grid_y_to_is_vdd_up_[y1] = 1;
  else
    grid_y_to_is_vdd_up_[y1] = 0;
}

void
Grid::markPowerMetalLayer(int layer_index, const LGRect* rect)
{
  auto [x1, y1, x2, y2] = getEndPointIndexExpanded(rect);

  // printf("(%d, %d) - (%d, %d): Layer: %d\n", x1, y1, x2, y2, layer_index);

  for(int i = x1; i <= x2; i++)
  {
    for(int j = y1; j <= y2; j++)
    {
      auto& pixel = pixels_[i][j];
      if(pixel.power_metal_layer == -1)
        pixel.power_metal_layer = layer_index;
      else
        pixel.power_metal_layer = std::min(layer_index, pixel.power_metal_layer);
      // TODO: What if multiple metals cover this pixel?
    }
  }
}

std::pair<PixelStatus, std::shared_ptr<LGRegion>> 
Grid::checkPixelStatus(const LGCell& cell) const
{
  const auto [x1, y1, x2, y2] = getEndPointIndexExpanded(&cell);

  for(int i = x1; i <= x2; i++)
  {
    for(int j = y1; j <= y2; j++)
    {
      const auto& pixel = pixels_[i][j];
      if(pixel.valid == false)
        return {PixelStatus::BLOCKAGE, pixel.region};
      else if(pixel.group_index != cell.getGroupIndex())
        return {PixelStatus::FENCE, pixel.region};
    }
  }

  return {PixelStatus::VALID, nullptr};
}

bool
Grid::isLegalMove(const LGCell& cell, int new_lx, int new_ly) const
{
  LGRect temp_rect(new_lx, new_ly, cell.getWidth(), cell.getHeight());

  const auto [x1, y1, x2, y2] = getEndPointIndexExpanded(&temp_rect);
  if(x1 < 0 or x2 >= x_size_ or y1 < 0 or y2 >= y_size_)
    return false;

  for(int i = x1; i <= x2; i++)
  {
    for(int j = y1; j <= y2; j++)
    {
      const auto& pixel = pixels_[i][j];
      if(pixel.valid == false)
        return false;
      if(pixel.group_index != cell.getGroupIndex())
        return false;
    }
  }

  return true;
}

}
