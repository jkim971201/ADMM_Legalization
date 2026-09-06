#ifndef GRID_H
#define GRID_H

#include <set>
#include <optional>
#include <vector>
#include <array>
#include <unordered_map>

#include "legalizer/Legalizer.h"

#include "db/dbSite.h"
#include "db/dbGroup.h"
#include "db/dbDatabase.h"

#include "objects/LGRegion.h"
#include "objects/LGCell.h"

namespace db {
  class dbRow;
  class dbDatabase;
}

namespace legalizer
{

enum class PixelStatus
{
  VALID,
  BLOCKAGE,
  FENCE,
  OUT_OF_CORE,
  POWER_RAIL_MISALIGN
};

struct Pixel
{
  bool valid;

  int h_overlap; // horizontal overlap with fence region
  int v_overlap; // vertical   overlap with fence region

  int power_metal_layer;
  int group_index;
  std::shared_ptr<LGRegion> region;

  Pixel() : valid(false), power_metal_layer(-1), 
            group_index(k_default_group_index), h_overlap(0), v_overlap(0) {}
};

class Grid : public LGRect
{
  public:

    Grid(std::shared_ptr<db::dbDatabase> db);
    
    int getSiteWidth() const { return site_width_; }
    int getRowHeight() const { return row_height_; }

    bool isLegalMove(const LGCell& cell, int new_lx, int new_ly) const;

    void markGroupRegion(const std::shared_ptr<LGRegion> region);
    void markNonGroupRegion(const std::shared_ptr<LGRegion> region);
    void markBlockageRegion(const std::shared_ptr<LGRegion> region);
    void markPowerMetalLayer(int layer_index, const LGRect* rect);
    void markDbRow(const db::dbRow* db_row);

    std::pair<PixelStatus, std::shared_ptr<LGRegion>> checkPixelStatus(const LGCell& cell) const;

    std::optional<std::pair<int, int>> findBestLegalPointNearby(
      const LGCell& cell, int real_x, int real_y) const;

    int getXSize() const { return x_size_; }
    int getYSize() const { return y_size_; }

    const std::vector<int>& getGridYToIsVddUp() const { return grid_y_to_is_vdd_up_; }

    const std::vector<std::vector<Pixel>>& getPixels() const { return pixels_; }
          std::vector<std::vector<Pixel>>& getPixels()       { return pixels_; }

    int getGridStartX(int dbu_x) const;
    int getGridStartY(int dbu_y) const;

    int getGridEndX(int dbu_x) const;
    int getGridEndY(int dbu_y) const;

    int getDbuX(int grid_x) const;
    int getDbuY(int grid_y) const;

    std::array<int, 4> getEndPointIndexExpanded(const LGRect* rect) const;
    std::array<int, 4> getEndPointIndexShrunk(const LGRect* rect) const;

  private:

    void initNonGroupGrid(std::shared_ptr<db::dbDatabase> db);

    int site_width_;
    int row_height_;

    int x_size_;
    int y_size_;

    std::vector<int> grid_y_to_is_vdd_up_;
    std::vector<std::vector<Pixel>> pixels_;
};

}

#endif
