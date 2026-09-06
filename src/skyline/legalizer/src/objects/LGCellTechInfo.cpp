#include "LGCellTechInfo.h"
#include "db/dbMacro.h"

#include <limits>

namespace legalizer
{

LGCellTechInfo::LGCellTechInfo(int id, const db::dbMacro* db_macro)
  : id_(id), ledge_type_(0), redge_type_(0), has_ground_at_bottom_(true), db_macro_(db_macro)
{
  ledge_type_ = db_macro->getLEdgeType();
  redge_type_ = db_macro->getREdgeType();

  int bottom_vss_y = std::numeric_limits<int>::max();
  int bottom_vdd_y = std::numeric_limits<int>::max();
  const auto db_mterms = db_macro->getMTerms();
  for(const auto& mterm : db_mterms)
  {
    const auto& ports_this_mterm = mterm->ports();
    for(const auto& port : ports_this_mterm)
    {
      const std::vector<std::pair<int, int>>& shape = port->getShape();
      for(const auto& [offsetX_dbu, offsetY_dbu] : shape)
      {
        if(mterm->usage() == PinUsage::POWER)
          bottom_vdd_y = std::min(bottom_vdd_y, offsetY_dbu);
        else if(mterm->usage() == PinUsage::GROUND)
          bottom_vss_y = std::min(bottom_vss_y, offsetY_dbu);
      }
    }
  }

  if(bottom_vss_y > bottom_vdd_y)
    has_ground_at_bottom_ = false;
  else
    has_ground_at_bottom_ = true;
}

std::string_view
LGCellTechInfo::getName() const 
{ 
  return db_macro_->name();
}

}
