#ifndef LG_CELL_TECH_INFO_H
#define LG_CELL_TECH_INFO_H

#include <map>

#include "LGRect.h"

namespace db {
  class dbMacro;
}

namespace legalizer
{

class LGCellTechInfo
{
  public:

    LGCellTechInfo(int id, const db::dbMacro* db_macro);

    int getID() const { return id_; }

    int getLEdgeType() const { return ledge_type_; }
    int getREdgeType() const { return redge_type_; }

    bool hasGroundAtBottom() const { return has_ground_at_bottom_; }

    std::string_view getName() const;

    // Tuple of {LayerIndex - PinBBox}
    const std::vector<std::tuple<int, int, int>>& getPinShapes() const { return pin_shape_; }
          std::vector<std::tuple<int, int, int>>& getPinShapes()       { return pin_shape_; }

  private:

    int id_;
    int ledge_type_;
    int redge_type_;

    bool has_ground_at_bottom_;

    const dbMacro* db_macro_;
    std::vector<std::tuple<int, int, int>> pin_shape_;
};

}

#endif
