#include "GPPin.h"
#include "GPCell.h"

#include "db/dbMTerm.h"

namespace skyplace
{

GPPin::GPPin()
  : id_       (0), 
    cx_       (0), 
    cy_       (0),
    offsetX_  (0), 
    offsetY_  (0),
    isIO_     (false),
    isMinPinX_(false), 
    isMinPinY_(false),
    isMaxPinX_(false), 
    isMaxPinY_(false),
    net_      (nullptr), 
    cell_     (nullptr),
    bterm_    (nullptr),   
    iterm_    (nullptr)
{}

GPPin::GPPin(dbITerm* iterm, int id, float scaling_factor)
  : GPPin()
{
  id_      = id;
  iterm_   = iterm;
  offsetX_ = static_cast<float>(iterm->getMTerm()->cx());
  offsetY_ = static_cast<float>(iterm->getMTerm()->cy());

  offsetX_ /= scaling_factor;
  offsetY_ /= scaling_factor;
}

GPPin::GPPin(dbBTerm* bterm, int id, float scaling_factor) 
  : GPPin()
{
  id_    = id;
  bterm_ = bterm;
  isIO_  = true;
  float bterm_cx = static_cast<float>(bterm->cx());
  float bterm_cy = static_cast<float>(bterm->cy());
  float bterm_lx = static_cast<float>(bterm->lx());
  float bterm_ly = static_cast<float>(bterm->ly());

  bterm_cx /= scaling_factor;
  bterm_cy /= scaling_factor;
  bterm_lx /= scaling_factor;
  bterm_ly /= scaling_factor;

  cx_      = bterm_cx;
  cy_      = bterm_cy;
  offsetX_ = bterm_cx - bterm_lx;
  offsetY_ = bterm_cy - bterm_ly;
}

void
GPPin::updatePinLocation(GPCell* cell) 
{
  // if we use cell_ instead of given GPCell*,
  // we have to check whether cell_ is not nullptr 
  // everytime we call this function
  cx_ = cell->cx() + offsetX_;
  cy_ = cell->cy() + offsetY_;
}

}
