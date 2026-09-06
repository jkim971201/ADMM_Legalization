#ifndef GP_RECT_H
#define GP_RECT_H

namespace skyplace
{

class GPRect
{
  public:

    GPRect();
    GPRect(int lx, int ly, int ux, int uy);

    float lx() const;
    float ly() const;
    float ux() const;
    float uy() const;

    float dx() const;
    float dy() const;

    float cx() const;
    float cy() const;

    float area() const;

    void setLx(float lx);
    void setLy(float ly);
    void setUx(float ux);
    void setUy(float uy);

  protected:

    float lx_;
    float ly_;
    float ux_;
    float uy_;
};

}

#endif
