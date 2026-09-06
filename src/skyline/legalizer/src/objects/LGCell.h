#ifndef LG_CELL_H
#define LG_CELL_H

#include <string>
#include <vector>
#include <functional>
#include <map>
#include <memory>

#include "LGRect.h"

#include "db/dbTypes.h"
#include "db/dbInst.h"

namespace legalizer
{

using namespace db;

class LGCellTechInfo;

class LGCell : public LGRect
{
  public:

    // Constructor for Regular Standard Cells
    LGCell(dbInst* inst);
      
    // Getters
    dbInst* getDbInst() const { return dbinst_; }
    int getInitLx() const { return init_lx_; }
    int getInitLy() const { return init_ly_; }
    int getInitUx() const { return init_lx_ + w_; }
    int getInitUy() const { return init_ly_ + h_; }
    int getInitCx() const { return init_lx_ + w_ / 2; }
    int getInitCy() const { return init_ly_ + h_ / 2; }

    int getDisp() const;

    int getLEdgeType() const;
    int getREdgeType() const; 

    double getWeight() const { return weight_; }

    Orient getOrient() const { return orient_; }

    bool isVddUp() const { return orient_ == N || orient_ == FN; }

    bool hasGroup() const { return group_index_ != -1; }

    std::string_view getName() const { return dbinst_->name(); }

    int getGroupIndex() const { return group_index_; }

    const std::shared_ptr<LGCellTechInfo> getTechInfo() const { return tech_info_; }

    // Setters 
    void setWeight(double weight) { weight_ = weight; }

    void setOrient(Orient orient) { orient_ = orient; }

    void setGroupIndex(int id) { group_index_ = id; }

    void setTechInfo(std::shared_ptr<LGCellTechInfo> info) { tech_info_ = info; }

  private:

    // Initial Lx, Ly are constant.
    // These will be used to compute total displacement.
    // (e.g. Global Placement results)
    int init_lx_;
    int init_ly_;

    int group_index_;

    int m2pin_bgn_;
    int m2pin_end_;

    double weight_; 

    Orient orient_;

    dbInst* dbinst_;

    std::shared_ptr<LGCellTechInfo> tech_info_;
};

}

#endif
