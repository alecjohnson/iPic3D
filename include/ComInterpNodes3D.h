/***************************************************************************
  ComInterpNodes.h  -  Library to manage communication of field vales among processors
  -------------------
begin                : May 2008
copyright            : (C) KU Leuven
developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/

#ifndef ComInterpNodes_H
#define ComInterpNodes_H

//#include "ComBasic3D.h"
#include "ipic_fwd.h"

/** communicate ghost cells and sum the contribution with a index indicating the number of species*/
void communicateInterp(int nx, int ny, int nz, double ***vector, const VirtualTopology3D * vct);

#endif
