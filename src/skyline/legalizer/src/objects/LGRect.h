#ifndef LG_RECT_H
#define LG_RECT_H

#include "db/dbInst.h"

namespace legalizer
{

using namespace db;

class LGRect
{
  public:

    LGRect();
    LGRect(const dbInst* inst);
    LGRect(int lx, int ly, int w, int h);

    void setLx(int lx);
    void setLy(int ly);

    void setWidth(int w);
    void setHeight(int h);

    int getLx() const;
    int getUx() const;

    int getLy() const;
    int getUy() const;

    int getWidth() const;
    int getHeight() const;

    int64_t getArea() const;

  protected:

    int lx_;
    int ly_;

    int w_;
    int h_;
};

}

#endif
