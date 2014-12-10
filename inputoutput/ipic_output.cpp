#include "ipic_output.h"

// *** end Grid3DCU_methods ***

// *** CAgetter_methods ***

arr3_double CAgetter::get_no_ghosts(const_arr3_double inarr)
{
  for (int i = 1; i < nxc-1; i++)
  for (int j = 1; j < nyc-1; j++)
  for (int k = 1; k < nzc-1; k++)
    arr[i-1][j-1][k-1]=inarr[i][j][k];
  return arr;
}
arr3_double CAgetter::get_N2C_no_ghosts(const_arr3_double inarr)
{
  array3_double tmp(nxc,nyc,nzc);
  grid->interpN2C(tmp, inarr);

  for (int i = 1; i < nxc-1; i++)
  for (int j = 1; j < nyc-1; j++)
  for (int k = 1; k < nzc-1; k++)
    arr[i-1][j-1][k-1]=inarr[i][j][k];
  return arr;
}
arr3_double CAgetter::get_N2C_no_ghosts(int is, const_arr4_double inarr)
{
  array4_double tmp(ns,nxc,nyc,nzc);
  grid->interpN2C(tmp, is, inarr);

  for (int i = 1; i < nxc-1; i++)
  for (int j = 1; j < nyc-1; j++)
  for (int k = 1; k < nzc-1; k++)
    arr[i-1][j-1][k-1]=tmp[is][i][j][k];
  return arr;
}

// *** end CAgetter_methods ***
