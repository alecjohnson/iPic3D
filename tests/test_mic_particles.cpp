#define MICVEC_DEFINE_OUTPUT_OPERATORS
#include <iostream>
#include <omp.h>
#include <stdio.h>
//#include "time.h" // for clock_gettime()
#include <stdint.h> // for uint64_t
#include <stdlib.h> // rand()
#include <sys/time.h>
#include <micvec.h>
#include <assert.h>
#include "Alloc.h"
#include "mic_particles.h"
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"

using namespace iPic3D;

const int nxc = 10;
const int nyc = 12;
const int nzc = 6;
const int nxn = nxc+1;
const int nyn = nyc+1;
const int nzn = nzc+1;
const int DFIELD_3or4 = 4;
array4<pfloat> fieldForPcls(nxn, nyn, nzn, 2*DFIELD_3or4);

void init_fieldForPcls()
{
  const int nxn = fieldForPcls.dim1();
  const int nyn = fieldForPcls.dim2();
  const int nzn = fieldForPcls.dim3();
  const int nnn = fieldForPcls.dim4();
  dprint(fieldForPcls.get_size());
  dprint(nxn);
  dprint(nyn);
  dprint(nzn);
  dprint(nnn);
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    fieldForPcls[i][j][k][0] = i+.2;
    fieldForPcls[i][j][k][1] = j+.3;
    fieldForPcls[i][j][k][2] = k+.4;
    fieldForPcls[i][j][k][0+DFIELD_3or4] = i+.5;
    fieldForPcls[i][j][k][1+DFIELD_3or4] = j+.6;
    fieldForPcls[i][j][k][2+DFIELD_3or4] = k+.7;
  }
}

void test_get_field_components_for_cell()
{
  init_fieldForPcls();
  F64vec8 field_components0[8];
  F64vec8 field_components1[8];
  I32vec16 cx(0,0,0,0,0,0,0,0,0,3,2,1,0,4,3,2);
  get_field_components_for_cell(
    field_components0,
    field_components1,
    fieldForPcls,
    cx);
  std::cout << "cx=" << cx << endl;
  for(int i=0;i<8;i++)
  {
    std::cout << "field_components0[" << i << "]=" << field_components0[i] << endl;
  }
  for(int i=0;i<8;i++)
  {
    std::cout << "field_components1[" << i << "]=" << field_components1[i] << endl;
  }
}

int main()
{
  printf("#\n");
  printf("# calling test_get_field_components_for_cell:\n");
  printf("#\n");
  test_get_field_components_for_cell();
}
