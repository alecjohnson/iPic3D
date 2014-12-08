/*******************************************************************************************
  CG.h  -  Conjugate Gradient Solver
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef CG_H
#define CG_H

typedef void (*CG_CALLBACK) (double *, double *, void **);

bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol,
  CG_CALLBACK FunctionImage, void ** registered_data);

#endif
