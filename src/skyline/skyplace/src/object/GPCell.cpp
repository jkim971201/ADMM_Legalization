#include "GPCell.h"
#include "GPPin.h"

namespace skyplace
{

GPCell::GPCell()
  : id_           (      0),
    cluster_id_   (      0),
    cx_           (      0), 
    cy_           (      0),
    dx_           (      0), 
    dy_           (      0),
    dDx_          (      0), 
    dDy_          (      0),
    densityScale_ (      1),
    isMacro_      (  false), 
    isFixed_      (  false),
    isFiller_     (  false),
    dbInst_       (nullptr)
{}

GPCell::GPCell(dbInst* inst) : GPCell()
{
  dbInst_  = inst;
  isFixed_ = inst->isFixed();
  isMacro_ = inst->isMacro();

  if(!isMacro_)
  {
    dx_ = static_cast<float>(inst->dx());
    dy_ = static_cast<float>(inst->dy());
    cx_ = static_cast<float>(inst->lx() + dx_ / 2);
    cy_ = static_cast<float>(inst->ly() + dy_ / 2);
  }
  else
  {
    dx_ = static_cast<float>(inst->dx() + inst->haloR() + inst->haloL());
    dy_ = static_cast<float>(inst->dy() + inst->haloT() + inst->haloB());
    cx_ = static_cast<float>(inst->lx() - inst->haloL() + dx_ / 2);
    cy_ = static_cast<float>(inst->ly() - inst->haloB() + dy_ / 2);
  }
}

GPCell::GPCell(dbInst* inst, float scaling_factor) : GPCell()
{
  dbInst_  = inst;
  isFixed_ = inst->isFixed();
  isMacro_ = inst->isMacro();

  if(!isMacro_)
  {
    dx_ = static_cast<float>(inst->dx());
    dy_ = static_cast<float>(inst->dy());
    cx_ = static_cast<float>(inst->lx() + dx_ / 2);
    cy_ = static_cast<float>(inst->ly() + dy_ / 2);
  }
  else
  {
    dx_ = static_cast<float>(inst->dx() + inst->haloR() + inst->haloL());
    dy_ = static_cast<float>(inst->dy() + inst->haloT() + inst->haloB());
    cx_ = static_cast<float>(inst->lx() - inst->haloL() + dx_ / 2);
    cy_ = static_cast<float>(inst->ly() - inst->haloB() + dy_ / 2);
  }

  dx_ /= scaling_factor;
  dy_ /= scaling_factor;
  cx_ /= scaling_factor;
  cy_ /= scaling_factor;
}

// Constructor for Filler GPCell
GPCell::GPCell(float cx, float cy, float dx, float dy) : GPCell()
{
  cx_ = cx;
  cy_ = cy;
  dx_ = dx;
  dy_ = dy;

  isFixed_  = false;
  isMacro_  = false;
  isFiller_ = true;
}

void 
GPCell::setCenterLocation(float newCx, float newCy)
{
  cx_ = newCx;
  cy_ = newCy;

  for(auto& pin : pins_)
    pin->updatePinLocation(this);
}

void 
GPCell::setDensitySize(float dWidth, float dHeight, float dScale)
{
  dDx_ = dWidth;
  dDy_ = dHeight;
  densityScale_ = dScale;
}

}
