
#include "mpi.h"
#include "ComNodes3D.h"
#include "ComParser3D.h"
#include "ComBasic3D.h"
#include "BcFields3D.h"
#include "VCtopology3D.h"
#include "TimeTasks.h"
#include "ipic_defs.h"
#include "Alloc.h"
#include "debug.h"
#include "parallel.h"

/** communicate ghost cells (FOR NODES) */
void communicateNode(int nx, int ny, int nz, arr3_double _vector,
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
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
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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
/** communicate ghost cells (FOR NODES) */
void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector = _vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;

  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
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
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, BCs, vct);

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
void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  const int BCs[6] = {
    bcFaceXrght, bcFaceXleft,
    bcFaceYrght, bcFaceYleft,
    bcFaceZrght, bcFaceZleft};
  communicateNodeBC(nx, ny, nz, _vector, BCs, vct);
}
/** communicate ghost cells (FOR NODES) with particles BC*/
void communicateNodeBC_P(int nx, int ny, int nz, arr3_double _vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
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
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

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

// PARTICLES
/** SPECIES: communicate ghost cells */
void communicateNode_P(int nx, int ny, int nz, double*** vector,
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
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
  makeNodeFace(nx, ny, nz, vector,
    ghostXrghtFace, ghostXleftFace,
    ghostYrghtFace, ghostYleftFace,
    ghostZrghtFace, ghostZleftFace);
  communicateGhostFace(nyr*nzr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtFace, ghostXleftFace);
  communicateGhostFace(nxr*nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostYrghtFace, ghostYleftFace);
  communicateGhostFace(nxr*nyr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostZrghtFace, ghostZleftFace);

  parseFace(nx, ny, nz, vector,
    ghostXrghtFace, ghostXleftFace,
    ghostYrghtFace, ghostYleftFace,
    ghostZrghtFace, ghostZleftFace);

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
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
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
  parseEdgeZ(nx, ny, nz, vector,
    ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge,
    ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector,
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge,
    ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector,
    ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge,
    ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



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
  parseCorner(nx, ny, nz, vector,
    &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner,
    &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner,
    &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner,
    &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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

// 
/** communicate ghost cells (FOR CENTERS) */
void communicateCenter(int nx, int ny, int nz, arr3_double _vector,
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector = _vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;

  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
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
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
void communicateCenterBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6], const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, BCs, vct);


  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
}
// particles
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
void communicateCenterBoxStencilBC_P(int nx, int ny, int nz,
  arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, BCs, vct);


  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
}

void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  const int BCs[6] = {
    bcFaceXrght, bcFaceXleft,
    bcFaceYrght, bcFaceYleft,
    bcFaceZrght, bcFaceZleft};
  BCface(nx, ny, nz, vector, BCs, vct);
  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
}
void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct)
{
  communicateNodeBoxStencilBC(nx, ny, nz, _vector,
    BCs[0], BCs[1], BCs[2], BCs[3], BCs[4], BCs[5], vct);
}

void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  const int BCs[6] = {
    bcFaceXrght, bcFaceXleft,
    bcFaceYrght, bcFaceYleft,
    bcFaceZrght, bcFaceZleft};
  communicateNodeBoxStencilBC_P(nx, ny, nz, _vector, BCs, vct);
}
void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6],
  const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, BCs, vct);
  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
}


// /////////// communication + BC ////////////////////////////
void communicateCenterBC(int nx, int ny, int nz, arr3_double vector_,
  const int BCs[6], const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=vector_.fetch_arr3();
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
  communicateGhostFace(nyr*nzr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtFace, ghostXleftFace);
  communicateGhostFace(nxr*nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostYrghtFace, ghostYleftFace);
  communicateGhostFace(nxr*nyr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector,
    ghostXrghtFace, ghostXleftFace,
    ghostYrghtFace, ghostYleftFace,
    ghostZrghtFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
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
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(),
    vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(),
    vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nxr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(),
    vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(),
    ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector,
    ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge,
    ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector,
    ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge,
    ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector,
    ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge,
    ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



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
  parseCorner(nx, ny, nz, vector,
    &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner,
    &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner,
    &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner,
    &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, BCs, vct);

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
void communicateCenterBC(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  const int BCs[6] = {
    bcFaceXrght, bcFaceXleft,
    bcFaceYrght, bcFaceYleft,
    bcFaceZrght, bcFaceZleft};
  communicateCenterBC(nx, ny, nz, vector, BCs, vct);
}
// /////////// communication + BC ////////////////////////////
void communicateCenterBC_P(int nx, int ny, int nz, arr3_double _vector,
  const int BCs[6], const VirtualTopology3D * vct)
{
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;

  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(nyr) * (nzr)];
  double *ghostXleftFace = new double[(nyr) * (nzr)];
  double *ghostYrghtFace = new double[(nxr) * (nzr)];
  double *ghostYleftFace = new double[(nxr) * (nzr)];
  double *ghostZrghtFace = new double[(nxr) * (nyr)];
  double *ghostZleftFace = new double[(nxr) * (nyr)];
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
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((nyr) * (nzr), vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nxr) * (nzr), vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nxr) * (nyr), vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(nyr, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nzr, vct->getCartesian_rank(), vct->getYrght(), vct->getYleft(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nxr, vct->getCartesian_rank(), vct->getZrght(), vct->getZleft(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXrght(), vct->getXleft(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, BCs, vct);

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
void communicateCenterBC_P(int nx, int ny, int nz, arr3_double vector,
  int bcFaceXrght, int bcFaceXleft,
  int bcFaceYrght, int bcFaceYleft,
  int bcFaceZrght, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  const int BCs[6] = {
    bcFaceXrght, bcFaceXleft,
    bcFaceYrght, bcFaceYleft,
    bcFaceZrght, bcFaceZleft};
  communicateCenterBC_P(nx, ny, nz, vector, BCs, vct);
}
