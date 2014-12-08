#ifndef GMRES_new2_H
#define GMRES_new2_H

#include "ipicfwd.h"

//typedef void (EMfields3D::*FIELD_IMAGE) (double *, double *);
typedef void (*GMRES_CALLBACK) (double *, double *, void**);

void GMRES(GMRES_CALLBACK FunctionImage, double *xkrylov, int xkrylovlen,
  const double *b, int m, int max_iter, double tol, void** registered_data);
void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);

#endif
