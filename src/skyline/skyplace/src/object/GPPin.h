#ifndef GP_PIN_H
#define GP_PIN_H

#include "db/dbITerm.h"
#include "db/dbBTerm.h"

namespace skyplace
{

using namespace db;

class GPNet;
class GPCell;

class GPPin
{
  public:

    GPPin();
    GPPin(dbBTerm* bterm, int id, float s); // Constructor for dbBTerm
    GPPin(dbITerm* iterm, int id, float s); // Constructor for dbITerm

    // Getters
    int     id() const { return id_;   }
    float   cx() const { return cx_;   }
    float   cy() const { return cy_;   }
    bool  isIO() const { return isIO_; }

    GPNet*   net() const { return net_;  }
    GPCell* cell() const { return cell_; } 
    // returns nullptr if GPPin is made from dbBTerm (== IO pin)

    float offsetX() const { return offsetX_; }
    float offsetY() const { return offsetY_; }

    bool isMinPinX() const { return isMinPinX_; }
    bool isMinPinY() const { return isMinPinY_; }
    bool isMaxPinX() const { return isMaxPinX_; }
    bool isMaxPinY() const { return isMaxPinY_; }

    dbBTerm* getDbBTerm() { return bterm_; }
    dbITerm* getDbITerm() { return iterm_; }

    std::string_view getName() const 
    { 
      return isIO_ ? bterm_->name() : iterm_->name(); 
    }

    // Setters
    void setId(int id)               { id_ = id;     }
    void setNet (GPNet*   net)       { net_  =  net; }
    void setCell(GPCell* cell)       { cell_ = cell; }
    void setOffset(float x, float y) { offsetX_ = x; offsetY_ = y; }

    void updatePinLocation(GPCell* cell);

    void   setMinPinX() { isMinPinX_ =  true; } 
    void   setMinPinY() { isMinPinY_ =  true; }
    void   setMaxPinX() { isMaxPinX_ =  true; }
    void   setMaxPinY() { isMaxPinY_ =  true; }
    void unsetMinPinX() { isMinPinX_ = false; }
    void unsetMinPinY() { isMinPinY_ = false; }
    void unsetMaxPinX() { isMaxPinX_ = false; }
    void unsetMaxPinY() { isMaxPinY_ = false; }

  private:

    int   id_;
    float cx_;
    float cy_;
    bool  isIO_;

    float offsetX_;
    float offsetY_;

    bool isMinPinX_;
    bool isMinPinY_;
    bool isMaxPinX_;
    bool isMaxPinY_;

    GPNet* net_;
    GPCell* cell_;

    dbBTerm* bterm_;
    dbITerm* iterm_;
};

}

#endif
