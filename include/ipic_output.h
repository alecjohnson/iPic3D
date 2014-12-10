#ifndef ipic_output_h
#define ipic_output_h
#include "Alloc.h"

// general methods used for output

// get cell-centered array by interpolating and stripping ghosts as needed
class CarrGetter
{
  private:
    array3_double arr;
    Grid3DCU*grid;

    CarrGetter(Grid3DCU*grid_):
      grid(grid_),
      arr(grid->get_nxc_r(),
          grid->get_nyc_r(),
          grid->get_nzc_r())
    {}
  public:
    arr3_double get_no_ghosts(const_arr3_double inarr);
    arr3_double get_N2C_no_ghosts(const_arr3_double inarr);
    arr3_double get_N2C_no_ghosts(int is, const_arr4_double inarr);
};

#endif // ipic_output_h
