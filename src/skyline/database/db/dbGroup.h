#ifndef DB_GROUP_H
#define DB_GROUP_H

#include <set>
#include <string>
#include <vector>

namespace db
{

class dbRegion;
class dbInst;

// FIXME
// At this moment, this can only handle 
// PatternName style group member (e.g. PatternName* ).
class dbGroup
{
  public:

    dbGroup(std::string_view name) : name_(name), pattern_style_(false) {}
  
    // Setters
    void setRegion(const dbRegion* region) { region_ = region; }
    void addPatternName(const std::string_view pattern) 
    { 
      pattern_style_ = true; 
      std::string new_pattern_name(pattern);
      pattern_set_.insert(new_pattern_name);
    }
  
    void addInst(dbInst* inst) { group_members_.push_back(inst); }

    // Getters
    const dbRegion* getRegion() const { return region_; }
    std::string_view getName() const { return name_; }
    bool isPatternStyle() const { return pattern_style_; }
    const std::set<std::string>& getPatternSet() const { return pattern_set_; }

          std::vector<dbInst*> getMembers()       { return group_members_; }
    const std::vector<dbInst*> getMembers() const { return group_members_; }

  private:
    
    std::string name_;

    bool pattern_style_;

    const dbRegion* region_;
    std::vector<dbInst*> group_members_; 

    std::set<std::string> pattern_set_;
};

}

#endif
