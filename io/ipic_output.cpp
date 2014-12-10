#include "ipic_output.h"
#include "Grid3DCU.h"

// *** end Grid3DCU_methods ***

// *** CarrGetter_methods ***

CarrGetter::CarrGetter(Grid3DCU*grid_, int ns_):
  grid(grid_),
  ns(ns_),
  nxc_r(grid->get_nxc_r()),
  nyc_r(grid->get_nyc_r()),
  nzc_r(grid->get_nzc_r()),
  nxc(grid->get_nxc()),
  nyc(grid->get_nyc()),
  nzc(grid->get_nzc()),
  arr(nxc_r, nyc_r, nzc_r)
{}

arr3_double CarrGetter::get_no_ghosts(const_arr3_double inarr)
{
  for (int i = 1; i <= nxc_r; i++)
  for (int j = 1; j <= nyc_r; j++)
  for (int k = 1; k <= nzc_r; k++)
    arr[i-1][j-1][k-1]=inarr[i][j][k];
  return arr;
}
arr3_double CarrGetter::get_N2C_no_ghosts(const_arr3_double inarr)
{
  array3_double tmp(nxc,nyc,nzc);
  grid->interpN2C(tmp, inarr);

  for (int i = 1; i <= nxc_r; i++)
  for (int j = 1; j <= nyc_r; j++)
  for (int k = 1; k <= nzc_r; k++)
    arr[i-1][j-1][k-1]=inarr[i][j][k];
  return arr;
}
arr3_double CarrGetter::get_N2C_no_ghosts(int is, const_arr4_double inarr)
{
  array4_double tmp(ns,nxc,nyc,nzc);
  grid->interpN2C(tmp, is, inarr);

  for (int i = 1; i <= nxc_r; i++)
  for (int j = 1; j <= nyc_r; j++)
  for (int k = 1; k <= nzc_r; k++)
    arr[i-1][j-1][k-1]=tmp[is][i][j][k];
  return arr;
}

// *** end CarrGetter_methods ***
