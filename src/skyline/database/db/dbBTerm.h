#ifndef DB_BTERM_H
#define DB_BTERM_H

#include <string>
#include <cstdlib>
#include <vector>

#include "dbTypes.h"
#include "dbLayer.h"
#include "dbRect.h"
#include "dbNet.h"

namespace db
{

class dbNet;
class dbBTerm;

class dbBTermPort : public dbRect
{
  public:

    // Constructor
    dbBTermPort()
      : offsetLx_ (0),
        offsetLy_ (0),
        offsetUx_ (0), 
        offsetUy_ (0),
        orient_   (Orient::N),
        status_   (Status::PLACED)
    {}

    void setLocation();

    // Setters
    void setOffsetLx(int offsetLx) { offsetLx_ = offsetLx; }
    void setOffsetLy(int offsetLy) { offsetLy_ = offsetLy; }
    void setOffsetUx(int offsetUx) { offsetUx_ = offsetUx; }
    void setOffsetUy(int offsetUy) { offsetUy_ = offsetUy; }

    void setOrient(Orient orient) { orient_ = orient; } 
    void setStatus(Status status) { status_ = status; } 

    void setOrigX(int val) { origX_ = val; }
    void setOrigY(int val) { origY_ = val; }

    void setBTerm(dbBTerm* bterm) { bterm_ = bterm; }

    // Getters
    int origX() const { return origX_; }
    int origY() const { return origY_; }

    int offsetLx() const { return offsetLx_; }
    int offsetLy() const { return offsetLy_; }
    int offsetUx() const { return offsetUx_; }
    int offsetUy() const { return offsetUy_; }

    bool isFixed()  const { return status_ == Status::FIXED;  }
    bool isCover()  const { return status_ == Status::COVER;  }
    bool isPlaced() const { return status_ == Status::PLACED; }

    Orient orient() const { return orient_;    }
    Status status() const { return status_;    }

    const dbBTerm* getBTerm() const { return bterm_; }
          dbBTerm* getBTerm()       { return bterm_; }

  private:

    int origX_;
    int origY_;

    int offsetLx_;
    int offsetLy_;
    int offsetUx_;
    int offsetUy_;

    Orient orient_;
    Status status_;

    dbBTerm* bterm_;
};

class dbBTerm
{
  public:

    dbBTerm();

    void print() const;

    // Setters
    void setName(const std::string& name) { name_ = name; }
    void setDirection (PinDirection dir) { direction_ = dir; }

    void setNet(dbNet* net) { net_ = net; }
    void addPort(dbBTermPort* port);

    // Getters
    const std::string& name() const { return name_; }
    PinDirection direction() const { return direction_; }

    const dbNet* net() const { return net_; }
    const std::vector<dbBTermPort*>& ports() const { return ports_; }
          std::vector<dbBTermPort*>& ports()       { return ports_; }

    // BBox
    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (ux_ + lx_) / 2; }
    int cy() const { return (uy_ + ly_) / 2; }
    int dx() const { return ux_ - lx_; }
    int dy() const { return uy_ - ly_; }

  private:

    int lx_; 
    int ly_;
    int ux_;
    int uy_;

    std::string name_;
    PinDirection direction_;

    dbNet* net_;
    std::vector<dbBTermPort*> ports_;
};

}

#endif
