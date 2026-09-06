#ifndef GP_BIN_H
#define GP_BIN_H

#include "GPRect.h"

namespace skyplace
{

class GPBin : public GPRect
{
  public:
    GPBin();
    GPBin(int row, int col, 
          float lx, float ly, 
          float ux, float uy, 
          float targetDensity);

    int row() const { return row_; }
    int col() const { return col_; }

    float density      () const { return density_;       }
    float targetDensity() const { return target_density_; }

    // Local Lagrange Multipler
    float lambda()        const { return lambda_; }

    float potential()     const { return electro_potential_; }

    float electroForceX() const { return electro_force_x_;    }
    float electroForceY() const { return electro_force_y_;    }

    void setLambda(float lambda) { lambda_ = lambda; }

    void setElectroPotential(float potential) { electro_potential_ = potential; }

    float fixedArea   () const { return fixed_area_;   }
    float movableArea () const { return movable_area_; }
    float fillerArea  () const { return filler_area_;  }

    void setDensity     (float density     ) { density_      = density;       }
    void setFixedArea   (float fixed_area  ) { fixed_area_   = fixed_area;    }
    void setMovableArea (float movable_area) { movable_area_ = movable_area;  }
    void setFillerArea  (float filler_area ) { filler_area_  = filler_area;   }

    void addFixedArea   (float fixed_area  ) { fixed_area_   += fixed_area;   }
    void addMovableArea (float movable_area) { movable_area_ += movable_area; }
    void addFillerArea  (float filler_area ) { filler_area_  += filler_area;  }

  private:

    int row_;
    int col_;

    float lambda_;

    float fixed_area_;
    float movable_area_;
    float filler_area_;

    float density_;
    float target_density_;

    float electro_potential_;
    float electro_force_x_;
    float electro_force_y_;
};

}

#endif
