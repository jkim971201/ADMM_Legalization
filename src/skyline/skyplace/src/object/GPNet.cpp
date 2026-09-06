#include <cstdio>
#include <numeric>

#include "GPNet.h"
#include "GPPin.h"

namespace skyplace
{

// GPNet //
GPNet::GPNet()
  : id_      (0),
    weight_  (1.0),
    dbNet_   (nullptr),
    GPRect   (0, 0, 0, 0)
{}

GPNet::GPNet(dbNet* net, int id) : GPNet()
{
  id_     =  id;
  weight_ = 1.0;
  dbNet_  = net;
}

void
GPNet::updateBBox()
{
  // To detect an error,
  // We initilize them as huge number.
  // so that an un-initilized net will
  // make total HPWL invalid.
  float new_lx = std::numeric_limits<float>::max();
  float new_ly = std::numeric_limits<float>::max();
  float new_ux = 0.0;
  float new_uy = 0.0;

  if(pins_.empty() == true)
  {
    printf("Warning - %s has no pins.\n", dbNetPtr()->name().c_str());
    new_lx = 0;
    new_ly = 0;
    new_ux = 0;
    new_uy = 0;
  }
  else
  {
    for(const auto& p : pins_)
    {
      new_lx = std::min(p->cx(), new_lx);
      new_ly = std::min(p->cy(), new_ly);
      new_ux = std::max(p->cx(), new_ux);
      new_uy = std::max(p->cy(), new_uy);
    }
  }
  
  setLx(new_lx);
  setLy(new_ly);
  setUx(new_ux);
  setUy(new_uy);
}

}
