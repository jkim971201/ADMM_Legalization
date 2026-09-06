#include "GPBin.h"

namespace skyplace
{

// GPBin //
GPBin::GPBin()
  : GPRect(0, 0, 0, 0),
    row_              (0), 
    col_              (0),
    lambda_           (1.0),
    fixed_area_       (0),
    movable_area_     (0),
    density_          (0),
    electro_force_x_  (0),
    electro_force_y_  (0),
    electro_potential_(0),
    target_density_   (0)
{}

GPBin::GPBin(
 int row, 
 int col, 
 float lx, 
 float ly, 
 float ux, 
 float uy,
 float target_density) 
  : GPRect(lx, ly, ux, uy),
    row_              (row), 
    col_              (col),
    lambda_           (1.0),
    fixed_area_       (0),
    movable_area_     (0),
    density_          (0),
    electro_force_x_  (0),
    electro_force_y_  (0),
    electro_potential_(0),
    target_density_   (target_density)
{}

}
