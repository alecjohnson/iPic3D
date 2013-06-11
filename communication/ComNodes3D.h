/***************************************************************************
  ComNodes.h  -  Library to manage communication of field values among processors
  -------------------
begin                : May 2008
copyright            : KUL Leuven
developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/

#ifndef ComNodes3D_H
#define ComNodes_H

#include "ComBasic3D.h"
#include "../utility/TimeTasks.h"
// boundary condition for fields
#include "../bc/BcFields3D.h"

/** communicate ghost cells (FOR NODES) */
inline void communicateNode(int nx, int ny, int nz, double ***vector, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
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
  timeTasks.addto_communicate();

}
/** communicate ghost cells (FOR NODES) */
inline void communicateNodeBC(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

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
  timeTasks.addto_communicate();

}
/** communicate ghost cells (FOR NODES) with particles BC*/
inline void communicateNodeBC_P(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
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
  timeTasks.addto_communicate();

}

/** SPECIES: communicate ghost cells */
inline void communicateNode(int nx, int ny, int nz, double ****vector, int ns, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);



/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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
  timeTasks.addto_communicate();

}                               // 

// PARTICLES
/** SPECIES: communicate ghost cells */
inline void communicateNode_P(int nx, int ny, int nz, double ****vector, int ns, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);



/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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
  timeTasks.addto_communicate();

}

// 
/** communicate ghost cells (FOR CENTERS) */
inline void communicateCenter(int nx, int ny, int nz, double ***vector, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);


/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
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
  timeTasks.addto_communicate();

}
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
inline void communicateCenterBoxStencilBC(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);


  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
  timeTasks.addto_communicate();
}
// particles
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
inline void communicateCenterBoxStencilBC_P(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);


  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
  timeTasks.addto_communicate();
}

// 


inline void communicateNodeBoxStencilBC(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
  timeTasks.addto_communicate();
}

inline void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
  // deallocate
  delete[]ghostXrghtFace;
  delete[]ghostXleftFace;
  delete[]ghostYrghtFace;
  delete[]ghostYleftFace;
  delete[]ghostZrghtFace;
  delete[]ghostZleftFace;
  timeTasks.addto_communicate();
}



/** SPECIES: communicate ghost cells */
inline void communicateCenter(int nx, int ny, int nz, double ****vector, int ns, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);


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
  timeTasks.addto_communicate();

}
// /////////// communication + BC ////////////////////////////
inline void communicateCenterBC(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  const int xFaceSize = (ny - 2) * (nz - 2);
  const int yFaceSize = (nx - 2) * (nz - 2);
  const int zFaceSize = (nx - 2) * (ny - 2);
  double *ghostXrghtFace = new double[xFaceSize*2];
  double *ghostXleftFace = new double[xFaceSize*2];
  double *ghostYrghtFace = new double[yFaceSize*2];
  double *ghostYleftFace = new double[yFaceSize*2];
  double *ghostZrghtFace = new double[zFaceSize*2];
  double *ghostZleftFace = new double[zFaceSize*2];
  // allocate 12 ghost cell Edges
  const int xEdgeSize = nx-2;
  const int yEdgeSize = ny-2;
  const int zEdgeSize = nz-2;
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[xEdgeSize*2];
  double *ghostXsameYrghtZleftEdge = new double[xEdgeSize*2];
  double *ghostXsameYleftZrghtEdge = new double[xEdgeSize*2];
  double *ghostXsameYrghtZrghtEdge = new double[xEdgeSize*2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[yEdgeSize*2];
  double *ghostXleftYsameZleftEdge = new double[yEdgeSize*2];
  double *ghostXrghtYsameZrghtEdge = new double[yEdgeSize*2];
  double *ghostXleftYsameZrghtEdge = new double[yEdgeSize*2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[zEdgeSize*2];
  double *ghostXrghtYrghtZsameEdge = new double[zEdgeSize*2];
  double *ghostXleftYleftZsameEdge = new double[zEdgeSize*2];
  double *ghostXleftYrghtZsameEdge = new double[zEdgeSize*2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner[2];
  double ghostXleftYrghtZrghtCorner[2];
  double ghostXrghtYleftZrghtCorner[2];
  double ghostXleftYleftZrghtCorner[2];
  double ghostXrghtYrghtZleftCorner[2];
  double ghostXleftYrghtZleftCorner[2];
  double ghostXrghtYleftZleftCorner[2];
  double ghostXleftYleftZleftCorner[2];

  MPI_Request requests[52];
  for(int i=0;i<52;i++) requests[i] = MPI_REQUEST_NULL;
  int request_idx = 0;
  int request_start = 0;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace(xFaceSize, 0, ghostXrghtFace, ghostXleftFace, &requests[4*request_idx++], vct);
  communicateGhostNowait(yFaceSize, 1, ghostYrghtFace, ghostYleftFace, &requests[4*request_idx++], vct);
  communicateGhostNowait(zFaceSize, 2, ghostZrghtFace, ghostZleftFace, &requests[4*request_idx++], vct);
  waitForRequests(requests,&request_start,request_idx);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  communicateGhostFace(ny - 2, 0, ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge, &requests[4*request_idx++], vct);
  communicateGhostFace(ny - 2, 0, ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge, &requests[4*request_idx++], vct);
  // Y-DIRECTION: X -> Y
  communicateGhostFace(nz - 2, 1, ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge, &requests[4*request_idx++], vct);
  communicateGhostFace(nz - 2, 1, ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge, &requests[4*request_idx++], vct);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nx - 2, 2, ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge, &requests[4*request_idx++], vct);
  communicateGhostFace(nx - 2, 2, ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge, &requests[4*request_idx++], vct);
  // parse
  waitForRequests(requests,&request_start,request_idx);
  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner, ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, 0, ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, &requests[4*request_idx++], vct);
  communicateGhostFace(1, 0, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner, &requests[4*request_idx++], vct);
  communicateGhostFace(1, 0, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner, &requests[4*request_idx++], vct);
  communicateGhostFace(1, 0, ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, &requests[4*request_idx++], vct);
  // parse
  assert(request_idx==13);
  waitForRequests(requests,&request_start,request_idx);
  parseCorner(nx, ny, nz, vector, ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner, ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

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
  timeTasks.addto_communicate();

}
// /////////// communication + BC ////////////////////////////
inline void communicateCenterBC_P(int nx, int ny, int nz, double ***vector, int bcFaceXrght, int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght, int bcFaceZleft, VirtualTopology3D * vct) {

  timeTasks.start_communicate();
  // allocate 6 ghost cell Faces
  double *ghostXrghtFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrghtFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrghtFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrghtZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrghtEdge = new double[nx - 2];
  double *ghostXsameYrghtZrghtEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrghtYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrghtYsameZrghtEdge = new double[ny - 2];
  double *ghostXleftYsameZrghtEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrghtYleftZsameEdge = new double[nz - 2];
  double *ghostXrghtYrghtZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrghtZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrghtYrghtZrghtCorner, ghostXleftYrghtZrghtCorner, ghostXrghtYleftZrghtCorner, ghostXleftYleftZrghtCorner;
  double ghostXrghtYrghtZleftCorner, ghostXleftYrghtZleftCorner, ghostXrghtYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrghtFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrghtFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrghtFace, ghostXleftFace, ghostYrghtFace, ghostYleftFace, ghostZrghtFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrghtFace, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrghtFace, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrghtFace, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYsameZrghtEdge, ghostXleftYsameZrghtEdge);
  // Y-DIRECTION: X -> Y
  MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrghtZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrghtYrghtZsameEdge, ghostXrghtYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrghtEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrghtZrghtEdge, ghostXsameYrghtZleftEdge);
  // parse
  MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrghtYrghtZsameEdge, ghostXleftYleftZsameEdge, ghostXrghtYleftZsameEdge, ghostXleftYrghtZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrghtYsameZrghtEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrghtEdge, ghostXrghtYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrghtZrghtEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrghtEdge, ghostXsameYrghtZleftEdge, &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner, &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner, &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner, &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZrghtCorner, &ghostXleftYrghtZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZrghtCorner, &ghostXleftYleftZrghtCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrghtYrghtZleftCorner, &ghostXleftYrghtZleftCorner);
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
  timeTasks.addto_communicate();

}

#endif
