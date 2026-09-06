#ifndef DB_REGION_H
#define DB_REGION_H

#include <string>
#include <vector>

#include "dbBox.h"
#include "dbTypes.h"

namespace db
{

class dbGroup;

class dbRegion 
{
  public:

    dbRegion(std::string_view name)
    {
      region_type_ = RegionType::FENCE;
      group_ = nullptr;
      name_ = name;
    }
  
    // Setters
    void setRegionType(RegionType type) { region_type_ = type; }
    void addRegionBox(int lx, int ly, int ux, int uy)
    {
      region_boxes_.push_back(dbBox(lx, ly, ux, uy));
    }
  
    void setGroup(dbGroup* group) { group_ = group; }

    // Getters
    RegionType getRegionType() const { return region_type_; }
    std::string_view getName() const { return name_; }
    const std::vector<dbBox>& getRegionBoxes() const { return region_boxes_; }

          dbGroup* getGroup()       { return group_; }
    const dbGroup* getGroup() const { return group_; }

  private:
    
    std::string name_;
    RegionType region_type_;
    std::vector<dbBox> region_boxes_;
    dbGroup* group_;
};

}

#endif
