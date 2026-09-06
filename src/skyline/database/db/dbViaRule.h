#ifndef DB_VIARULE_H
#define DB_VIARULE_H

#include <string>
#include <string_view>
#include <vector>

namespace db
{

class dbViaRule
{
  public:

    dbViaRule(const char* name) : name_(std::string(name)) {}

    // Setters
    
    // Getters
    std::string_view getName() const { return name_; }

  private:

		std::string name_;

		// NOTE
		// Need to make dbViaRuleLayerRule (naming?)
};

}

#endif
