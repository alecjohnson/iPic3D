/*******************************************************************************************
  CG.h  -  Conjugate Gradient Solver
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef CG_H
#define CG_H

#include "ipicfwd.h"

// These declarations are currently needed because Field is not anymore a class
// CG needs a pointer to the function that solves the fields.
// 
// To avoid changing all the code we typedef Field as of type EMfields3D (which is
// not derived anymore from Field). This will be improved in future releases.

class EMfields3D;
typedef EMfields3D Field;
typedef void (Field::*FIELD_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);
typedef void (*GENERIC_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);

bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, FIELD_IMAGE FunctionImage, Grid * grid, VirtualTopology3D * vct, Field * field);
bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, GENERIC_IMAGE FunctionImage, Grid * grid, VirtualTopology3D * vct);

#endif
