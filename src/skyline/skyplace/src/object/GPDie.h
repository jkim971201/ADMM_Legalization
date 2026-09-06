#ifndef GP_DIE_H
#define GP_DIE_H

#include "GPRect.h"

namespace skyplace
{

// In the current implementation,
// lx ly ux uy are coodinates of Core,
// where we just assume that die ly ly is (0, 0).
// So, if (lx, ly) of DEF DIEAREA is not (0, 0),
// then this will make bug.
class GPDie : public GPRect
{
  public:

    GPDie(int _lx, int _ly, int _ux, int _uy) 
      : GPRect(_lx, _ly, _ux, _uy)
    {}

  private:
};

}

#endif
