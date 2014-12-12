
#ifndef __PARALLELIO_H__
#define __PARALLELIO_H__

#ifdef USEH5HUT
#  include "H5hut-io.h"
#endif

#ifdef PHDF5
#  include "phdf5.h"
#endif

#include "ipic_fwd.h"

void WriteFieldsH5hut(int nspec,
  const CollectiveIO *col,
  const VCtopology3D *vct,
  const Grid3DCU *grid,
  const EMfields3D *EMf,
  int cycle);
void WritePartclH5hut(int nspec,
  const CollectiveIO *col,
  const VCtopology3D *vct,
  const Grid3DCU *grid,
  const Particles3Dcomm *part,
  int cycle);

void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid);
void ReadFieldsH5hut(int nspec, EMfields3D *EMf,       Collective *col, VCtopology3D *vct, Grid3DCU *grid);

void WriteOutputParallel(
  const CollectiveIO *col,
  const VCtopology3D *vct,
  const Grid3DCU *grid,
  const Particles3Dcomm *part,
  const SpeciesMoms *speciesMoms,
  const EMfields3D *EMf,
  int cycle);

#endif
