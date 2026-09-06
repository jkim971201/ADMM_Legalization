#ifndef GP_CELL_H
#define GP_CELL_H

#include <vector>

#include "db/dbInst.h"

namespace skyplace
{

using namespace db;

class GPPin;

class GPCell 
{
  public:

    GPCell();
    GPCell(dbInst* inst);
    GPCell(dbInst* inst, float scaling_factor);
    GPCell(float cx, float cy, float dx, float dy);
    // Contructor for Filler Cell

    // Getters
    int   id() const { return id_;           }
    float lx() const { return cx_ - dx_ / 2; }
    float ly() const { return cy_ - dy_ / 2; }
    float ux() const { return cx_ + dx_ / 2; }
    float uy() const { return cy_ + dy_ / 2; }

    float cx() const { return cx_; }
    float cy() const { return cy_; }
    float dx() const { return dx_; }
    float dy() const { return dy_; }

    float dLx() const { return cx_ - dDx_ / 2; }
    float dLy() const { return cy_ - dDy_ / 2; }
    float dUx() const { return cx_ + dDx_ / 2; }
    float dUy() const { return cy_ + dDy_ / 2; }
    float dDx() const { return dDx_; }
    float dDy() const { return dDy_; }

    float area() const { return dx_ * dy_; }
    float densityScale() const { return densityScale_; }

    bool isFixed  () const { return isFixed_ ;   }
    bool isMacro  () const { return isMacro_ ;   }
    bool isFiller () const { return isFiller_;   }

    int  clusterID() const { return cluster_id_; }

    std::string_view getName() const { return dbInst_->name(); }

    const std::vector<GPPin*>& pins() const { return pins_; }

    // These will return nullptr,
    // if there is no corresponding dbInst (e.g. fillerCell)
    dbInst* dbInstPtr() const { return dbInst_;  }

    // Setters
    void setID       (int id) { id_  = id; }
    void setClusterID(int id) { cluster_id_ = id; }

    void setCenterLocation (float newCx,  float newCy);
    void setDensitySize    (float dWidth, float dHeight, float dScale);
    void setDensityScale   (float dScale) { densityScale_ = dScale; }

    void addNewPin(GPPin* pin) { pins_.push_back(pin); }

  private:

    dbInst*  dbInst_;

    // ID is required to compute Laplacian
    // Fixed and Movables seperately
    int id_;
    int cluster_id_;

    float cx_;
    float cy_;

    float dx_;
    float dy_;

    float dDx_;
    float dDy_;

    float densityScale_;

    bool isMacro_;
    bool isFixed_;
    bool isFiller_;

    std::vector<GPPin*> pins_;
};

}

#endif
