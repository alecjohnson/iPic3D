/***************************************************************************
  BcFields3D.h  -  Library to manage boundary conditions for fields 
  -------------------
begin                : Fri Jan 2009
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef BcFields_H
#define BcFields_H

#include "VCtopology3D.h"

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, double ***vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);
inline void BCface(int nx, int ny, int nz, double ***vector,
  const int BCs[6], const VirtualTopology3D * vct)
{
  BCface(nx, ny, nz, vector,
    BCs[0], BCs[1], BCs[2], BCs[3], BCs[4], BCs[5], vct);
}

/** set the boundary condition on boundaries */
void BCface_P(int nx, int ny, int nz, double ***vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct);
inline void BCface_P(int nx, int ny, int nz, double ***vector,
  const int BCs[6], const VirtualTopology3D * vct)
{
  BCface_P(nx, ny, nz, vector,
    BCs[0], BCs[1], BCs[2], BCs[3], BCs[4], BCs[5], vct);
}


/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, int ns, double ****vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct);

/** set the boundary condition on boundaries Particles*/
void BCface_P(int nx, int ny, int nz, int ns, double ****vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct);

#endif
