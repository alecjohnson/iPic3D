/***************************************************************************
  ComNodes.h  -  Library to manage communication of field values among processors
  -------------------
begin                : May 2008
copyright            : KUL Leuven
developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/

#ifndef ComNodes3D_H
#define ComNodes3D_H

#include "arraysfwd.h"
#include "ipic_fwd.h"

class Setting;

// boundary condition for fields
//#include "BcFields3D.h"

/** communicate ghost cells (FOR NODES) */
void communicateNode(int nx, int ny, int nz, arr3_double vector,
  const VirtualTopology3D * vct);

/** communicate ghost cells (FOR NODES) */
void communicateNodeBC(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
// version with fewer arguments
void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct);

/** communicate ghost cells (FOR NODES) with particles BC*/
void communicateNodeBC_P(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);

// PARTICLES
/** SPECIES: communicate ghost cells */
void communicateNode_P(int nx, int ny, int nz, double*** vector, const VirtualTopology3D * vct);

// 
/** communicate ghost cells (FOR CENTERS) */
void communicateCenter(int nx, int ny, int nz, arr3_double vector, const VirtualTopology3D * vct);

/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
//void communicateCenterBoxStencilBC(int nx, int ny, int nz, arr3_double vector,
//int bcFaceXright, int bcFaceXleft,
//int bcFaceYright, int bcFaceYleft,
//int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
void communicateCenterBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6], const VirtualTopology3D * vct);

// particles
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
void communicateCenterBoxStencilBC_P(int nx, int ny, int nz, arr3_double vector,
  const int BCs[6], const VirtualTopology3D * vct);

void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
// version with fewer arguments
void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct);

void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
// version with fewer arguments
void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double vector,
  const int BCs[6],
  const VirtualTopology3D * vct);

// /////////// communication + BC ////////////////////////////
void communicateCenterBC(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
// version with fewer arguments
void communicateCenterBC(int nx, int ny, int nz, arr3_double vector_,
  const int BCs[6], const VirtualTopology3D * vct);

// /////////// communication + BC ////////////////////////////
void communicateCenterBC_P(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
// version with fewer arguments
void communicateCenterBC_P(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6], const VirtualTopology3D * vct);

#endif
