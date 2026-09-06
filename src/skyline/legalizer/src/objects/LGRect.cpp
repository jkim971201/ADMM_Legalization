#include "LGRect.h"

namespace legalizer
{

LGRect::LGRect() : lx_(0), ly_(0), w_(0), h_(0) {}

LGRect::LGRect(const dbInst* inst)
{
  lx_ = inst->lx();
  ly_ = inst->ly();
  w_  = inst->dx();
  h_  = inst->dy();
}

LGRect::LGRect(int lx, int ly, int w, int h) : lx_(lx), ly_(ly), w_(w), h_(h) {}

void LGRect::setLx(int lx) { lx_ = lx; }
void LGRect::setLy(int ly) { ly_ = ly; }
void LGRect::setWidth(int w) { w_ = w; }
void LGRect::setHeight(int h) { h_ = h; }

int LGRect::getLx() const { return lx_; }
int LGRect::getUx() const { return lx_ + w_; }
int LGRect::getLy() const { return ly_; }
int LGRect::getUy() const { return ly_ + h_; }
int LGRect::getWidth() const { return w_; }
int LGRect::getHeight() const { return h_; }

int64_t LGRect::getArea() const { return int64_t(w_) * int64_t(h_); }

}
