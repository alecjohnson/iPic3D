#ifndef ipic_output_h
#define ipic_output_h
class Grid3DCU;
#include "Alloc.h"

// general methods used for output

// get cell-centered array by interpolating and stripping ghosts as needed
class CarrGetter
{
  private:
    array3_double arr;
    Grid3DCU*grid;
    const int ns;
    const int nxc_r;
    const int nyc_r;
    const int nzc_r;
    const int nxc;
    const int nyc;
    const int nzc;

  public:
    CarrGetter(Grid3DCU*grid_, int ns_);
    arr3_double get_no_ghosts(const_arr3_double inarr);
    arr3_double get_N2C_no_ghosts(const_arr3_double inarr);
    arr3_double get_N2C_no_ghosts(int is, const_arr4_double inarr);
};

#endif // ipic_output_h
