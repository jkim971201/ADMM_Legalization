#ifndef LG_REGION_H
#define LG_REGION_H

#include "LGRect.h"
#include <vector>

namespace legalizer
{

enum class RegionType
{
  DEFAULT,
  FENCE,
  BLOCKAGE
};

enum class RegionShape
{
  NORMAL = 0,
  THIN_HORIZONTAL = 1,
  THIN_VERTICAL = 2
};

class LGRegion : public LGRect
{
  public:

    // lx ly ux uy are in dbu 
    LGRegion(int id, int lx, int ly, int ux, int uy, int group_index, RegionType type)
      : id_(id), group_index_(group_index), type_(type), shape_(RegionShape::NORMAL), LGRect(lx, ly, ux - lx, uy - ly) {}

    int getIndex() const { return id_; }

    int getGroupIndex() const { return group_index_; }
    RegionType getType() const { return type_; }
    RegionShape getShape() const { return shape_; }

    void addAdjacentRegion(std::shared_ptr<LGRegion> region) { adjacent_regions_.push_back(region); }

    void setShape(RegionShape new_shape) { shape_ = new_shape; }

    bool isIsolated() const { return adjacent_regions_.empty(); }

  private:

    int id_;
    int group_index_;
    RegionType type_;
    RegionShape shape_;

    std::vector<std::shared_ptr<LGRegion>> adjacent_regions_;
};

}

#endif
