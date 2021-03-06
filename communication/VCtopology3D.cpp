#include "mpi.h"
#include "Alloc.h"
#include "Collective.h"
#include "VCtopology3D.h"
#include <iostream>
#include "MPIdata.h"
#include "debug.h"

using std::cout;
using std::endl;

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
VCtopology3D::VCtopology3D(const Collective& col) {
  // *******************************************
  // *******************************************
  // change these values to change the topology
  XLEN = col.getXLEN();
  YLEN = col.getYLEN();
  ZLEN = col.getZLEN();
  nprocs = XLEN * YLEN * ZLEN;
  // here you have to set the topology
  PERIODICX = col.getPERIODICX();
  PERIODICY = col.getPERIODICY();
  PERIODICZ = col.getPERIODICZ();
  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;

  reorder = 1;

  dims[0] = XLEN;
  dims[1] = YLEN;
  dims[2] = ZLEN;

  periods[0] = PERIODICX;
  periods[1] = PERIODICY;
  periods[2] = PERIODICZ;

  cVERBOSE = false;             // communication verbose ?

}





/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  */
void VCtopology3D::setup_vctopology(MPI_Comm old_comm) {
  // create a matrix with ranks, and neighbours for fields
  MPI_Cart_create(old_comm, 3, dims, periods, reorder, &CART_COMM);
  // create a matrix with ranks, and neighbours for Particles
  //MPI_Cart_create(old_comm, 3, dims, periods_P, reorder, &CART_COMM_P);
  // Why not the following line instead of the previous?  Was
  // this written in anticipation that a different number of MPI
  // processes would be used for fields versus for particles?
  // But the code has not been consistently written this way...
  //MPI_Cart_create(CART_COMM, 3, dims, periods_P, 0, &CART_COMM_P);
  // field Communicator
  if (CART_COMM != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
    MPI_Cart_coords(CART_COMM, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM, XDIR, RIGHT, &xleft_neighbor, &xright_neighbor);
    MPI_Cart_shift(CART_COMM, YDIR, RIGHT, &yleft_neighbor, &yright_neighbor);
    MPI_Cart_shift(CART_COMM, ZDIR, RIGHT, &zleft_neighbor, &zright_neighbor);
  }
  else {
    // previous check that nprocs = XLEN*YLEN*ZLEN should prevent reaching this line.
    eprintf("A process is thrown away from the new topology for fields.");
  }
  // Particles Communicator
  //if (CART_COMM_P != MPI_COMM_NULL) {
  //  int pcl_coordinates[3];
  //  int pcl_cartesian_rank;
  //  MPI_Comm_rank(CART_COMM_P, &pcl_cartesian_rank);
  //  MPI_Cart_coords(CART_COMM_P, pcl_cartesian_rank, 3, pcl_coordinates);
  //  
  // This seems to be assumed elsewhere in the code.
  // We need to eliminate this assumption.
  // This becomes important if we want
  // a running program to have more than one
  // mesh and processor topology.
  // The MPI rank is for system-level code, e.g.
  // to identify the process from which debug is coming.
  // The cartesian rank is for application-level code.
  assert_eq(cartesian_rank, MPIdata::get_rank());
  //  // should agree
  //  assert_eq(cartesian_rank,pcl_cartesian_rank);
  //  for(int dim=0;dim<3;dim++)
  //  {
  //    assert_eq(coordinates[dim],pcl_coordinates[dim]);
  //  }
  //
  //  MPI_Cart_shift(CART_COMM_P, XDIR, RIGHT, &xleft_neighbor_P, &xright_neighbor_P);
  //  MPI_Cart_shift(CART_COMM_P, YDIR, RIGHT, &yleft_neighbor_P, &yright_neighbor_P);
  //  MPI_Cart_shift(CART_COMM_P, ZDIR, RIGHT, &zleft_neighbor_P, &zright_neighbor_P);
  //}
  //else {
  //  // previous check that nprocs = XLEN*YLEN*ZLEN should prevent reaching this line.
  //  eprintf("A process is thrown away from the new topology for Particles.");
  //}

  _isPeriodicXlower = PERIODICX && isXlower();
  _isPeriodicXupper = PERIODICX && isXupper();
  _isPeriodicYlower = PERIODICY && isYlower();
  _isPeriodicYupper = PERIODICY && isYupper();
  _isPeriodicZlower = PERIODICZ && isZlower();
  _isPeriodicZupper = PERIODICZ && isZupper();

  _noXleftNeighbor = (getXleft() == MPI_PROC_NULL);
  _noXrghtNeighbor = (getXrght() == MPI_PROC_NULL);
  _noYleftNeighbor = (getYleft() == MPI_PROC_NULL);
  _noYrghtNeighbor = (getYrght() == MPI_PROC_NULL);
  _noZleftNeighbor = (getZleft() == MPI_PROC_NULL);
  _noZrghtNeighbor = (getZrght() == MPI_PROC_NULL);

  _isBoundaryProcess = 
    _noXleftNeighbor ||
    _noXrghtNeighbor ||
    _noYleftNeighbor ||
    _noYrghtNeighbor ||
    _noZleftNeighbor ||
    _noZrghtNeighbor;
}
/** destructor */
VCtopology3D::~VCtopology3D() {

}
/** print topology info */
void VCtopology3D::Print() {
  cout << endl;
  cout << "Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Processors grid: " << XLEN << "x" << YLEN << "x" << ZLEN << endl;
  cout << "Periodicity X: " << periods[0] << endl;
  cout << "Periodicity Y: " << periods[1] << endl;
  cout << "Periodicity Z: " << periods[2] << endl;
  cout << endl;
}
/** print cartesian rank of neighbors and coordinate of process */
void VCtopology3D::PrintMapping() {
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << "; Z = " << coordinates[2] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = " << yright_neighbor << "; zLeft = " << zleft_neighbor << "; zRight = " << zright_neighbor << endl;
  cout << endl;
}
