#include "LGCell.h"
#include "LGCellTechInfo.h"

namespace legalizer
{

// Constructor for Regular Standard Cells
LGCell::LGCell(dbInst* inst)
  : dbinst_     (inst),
    init_lx_    (inst->lx()),
    init_ly_    (inst->ly()),
    weight_     (1.0),
    group_index_(-1),
    m2pin_bgn_  (-1),
    m2pin_end_  (-1),
    orient_     (inst->orient()),
    LGRect      (inst)
{}

int LGCell::getDisp() const 
{ 
  return std::abs(lx_ - init_lx_) + std::abs(ly_ - init_ly_); 
}

int 
LGCell::getLEdgeType() const 
{ 
  return tech_info_->getLEdgeType(); 
}

int 
LGCell::getREdgeType() const 
{ 
  return tech_info_->getREdgeType(); 
}

}
