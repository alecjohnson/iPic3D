#include <mpi.h>
#include "ComInterpNodes3D.h"
#include "ComParser3D.h"
#include "ComBasic3D.h"
#include "ipicdefs.h"
#include "Alloc.h"
#include "VCtopology3D.h"

/** communicate and sum shared ghost cells */
void communicateInterp(int nx, int ny, int nz, double*** vector, const VirtualTopology3D * vct)
{
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[nyr*nzr];
  double *ghostXleftFace = new double[nyr*nzr];
  double *ghostYrghtFace = new double[nxr*nzr];
  double *ghostYleftFace = new double[nxr*nzr];
  double *ghostZrghtFace = new double[nxr*nyr];
  double *ghostZleftFace = new double[nxr*nyr];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nxr];
  double *ghostXsameYrghtZleftEdge = new double[nxr];
  double *ghostXsameYleftZrghtEdge = new double[nxr];
  double *ghostXsameYrghtZrghtEdge = new double[nxr];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[nyr];
  double *ghostXleftYsameZleftEdge = new double[nyr];
  double *ghostXrghtYsameZrghtEdge = new double[nyr];
  double *ghostXleftYsameZrghtEdge = new double[nyr];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nzr];
  double *ghostXrghtYrghtZsameEdge = new double[nzr];
  double *ghostXleftYleftZsameEdge = new double[nzr];
  double *ghostXleftYrghtZsameEdge = new double[nzr];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner;
  double ghostXleftYrghtZrghtCorner;
  double ghostXrghtYleftZrghtCorner;
  double ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner;
  double ghostXleftYrghtZleftCorner;
  double ghostXrghtYleftZleftCorner;
  double ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector,
    ghostXrghtFace, ghostXleftFace,
    ghostYrghtFace, ghostYleftFace,
    ghostZrghtFace, ghostZleftFace);
  // communication
  communicateGhostFace(nyr*nzr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtFace, ghostXleftFace);
  communicateGhostFace(nxr*nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostYrghtFace, ghostYleftFace);
  communicateGhostFace(nxr*nyr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostZrghtFace, ghostZleftFace);

  addFace(nx, ny, nz, vector,
    ghostXrghtFace, ghostXleftFace,
    ghostYrghtFace, ghostYleftFace,
    ghostZrghtFace, ghostZleftFace, vct);

/** prepare ghost cell Edge Y for communication: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace,
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge,
    ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace,
    ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge,
    ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace,
    ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge,
    ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nxr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  addEdgeZ(nx, ny, nz, vector,
    ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge,
    ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge, vct);
  addEdgeY(nx, ny, nz, vector,
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge,
    ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge, vct);
  addEdgeX(nx, ny, nz, vector,
    ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge,
    ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, vct);



  makeNodeCorner(nx, ny, nz,
    ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge,
    ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge,
    &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner,
    &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner,
    &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner,
    &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  addCorner(nx, ny, nz, vector,
    &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner,
    &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner,
    &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner,
    &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner, vct);


  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrghtZleftEdge;
  delete[]ghostXsameYleftZrghtEdge;
  delete[]ghostXsameYrghtZrghtEdge;
  // Y EDGE
  delete[]ghostXrghtYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrghtYsameZrghtEdge;
  delete[]ghostXleftYsameZrghtEdge;
  // Z EDGE
  delete[]ghostXrghtYleftZsameEdge;
  delete[]ghostXrghtYrghtZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrghtZsameEdge;
}
