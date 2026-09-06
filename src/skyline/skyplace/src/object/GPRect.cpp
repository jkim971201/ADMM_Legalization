#include "GPRect.h"

namespace skyplace
{

GPRect::GPRect()
  : lx_(0), ly_(0), ux_(0), uy_(0)
{}

GPRect::GPRect(int _lx, int _ly, int _ux, int _uy)
{
  lx_ = static_cast<float>(_lx);
  ly_ = static_cast<float>(_ly);
  ux_ = static_cast<float>(_ux);
  uy_ = static_cast<float>(_uy);
}

float
GPRect::lx() const
{
  return lx_;
}

float
GPRect::ly() const
{
  return ly_;
}

float
GPRect::ux() const
{
  return ux_;
}

float
GPRect::uy() const
{
  return uy_;
}

float
GPRect::dx() const
{
  return ux_ - lx_;
}

float
GPRect::dy() const
{
  return uy_ - ly_;
}

float 
GPRect::cx() const 
{ 
  return (lx_ + ux_) / 2.0;
}

float 
GPRect::cy() const 
{ 
  return (ly_ + uy_) / 2.0; 
}

float
GPRect::area() const
{
  return dx() * dy();
}

void
GPRect::setLx(float lx)
{
  lx_ = lx;
}

void
GPRect::setLy(float ly)
{
  ly_ = ly;
}

void
GPRect::setUx(float ux)
{
  ux_ = ux;
}

void
GPRect::setUy(float uy)
{
  uy_ = uy;
}

}
