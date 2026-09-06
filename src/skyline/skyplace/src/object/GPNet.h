#ifndef GP_NET_H
#define GP_NET_H

#include <vector>

#include "db/dbNet.h"

#include "GPRect.h"

namespace skyplace
{

using namespace db;

class GPPin;

// Inherit GPRect to describe
// the BBox of this net.
class GPNet : public GPRect
{
  public:

    GPNet();
    GPNet(dbNet* net, int id); 

    // Getters
    int   id()     const { return id_;          }
    int   deg()    const { return pins_.size(); }
    float hpwl()   const { return dx() + dy();  }
    float weight() const { return weight_;      }

    const std::vector<GPPin*>& pins() const { return pins_; }

    dbNet* dbNetPtr() const { return dbNet_; }

    // Setters
    void setId      (int id)       { id_ = id; }
    void setDbNet   (dbNet* dbnet) { dbNet_ = dbnet; }
    void setWeight  (float weight) { weight_ = weight; }
    void addNewPin  (GPPin* pin)     { pins_.push_back(pin); }

    void updateBBox();

  private:

    dbNet* dbNet_;

    int id_;
    float weight_;

    std::vector<GPPin*> pins_;
};

}

#endif
