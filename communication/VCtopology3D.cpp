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
  dims[0] = col.getXLEN();
  dims[1] = col.getYLEN();
  dims[2] = col.getZLEN();
  nprocs = dims[0] * dims[1] * dims[2];
  // here you have to set the topology
  periods[0] = col.getPERIODICX();
  periods[1] = col.getPERIODICY();
  periods[2] = col.getPERIODICZ();
  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;

  reorder = 1;

  cVERBOSE = false;             // communication verbose ?

}

// problem: MPI standard does not specify whether
// row-major or column-major indexing is used, as
// far as I know.
//void VCtopology::compute_rank_for_cart_coords(int i, int j, int k)
//{
//  return (i*dims[1]+j)*dims[2]+k;
//}
//void VCtopology::compute_cart_coords_for_rank(int coords[3], int r)
//{
//  coords[2] = r%dims[2];
//  r/=dims[2];
//  coords[1] = r%dims[1];
//  r/=dims[1];
//  coords[0] = r;
//}


/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  */
void VCtopology3D::setup_vctopology(MPI_Comm old_comm) {

  if(getNprocs() != MPIdata::get_nprocs())
  {
    if(MPIdata::get_rank() == 0)
    {
      cout << "Error: " << MPIdata::get_nprocs()
        << " processes can't be mapped to a "
        << getXLEN() << "x"
        << getYLEN() << "x"
        << getZLEN()
        << " matrix: Change XLEN, YLEN, or ZLEN in .inp file"
        << endl;
      MPIdata::instance().finalize_mpi();
      abort();
    }
  }
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
  if (CART_COMM != MPI_COMM_NULL)
  {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
    MPI_Cart_coords(CART_COMM, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM, XDIR, RIGHT,&_leftNeighbor[0],&_rghtNeighbor[0]);
    MPI_Cart_shift(CART_COMM, YDIR, RIGHT,&_leftNeighbor[1],&_rghtNeighbor[1]);
    MPI_Cart_shift(CART_COMM, ZDIR, RIGHT,&_leftNeighbor[2],&_rghtNeighbor[2]);

    // initialize ranks of diagonal neighbors
    //
    for(int i0=-1;i0<=1;i0++)
    for(int i1=-1;i1<=1;i1++)
    for(int i2=-1;i2<=1;i2++)
    {
      int coords[3] = {
        coordinates[0]+i0,
        coordinates[1]+i1,
        coordinates[2]+i2 };
      int rank;
      if((!periods[0] && (coords[0]<0 || coords[0]>=dims[0]))
       ||(!periods[1] && (coords[1]<0 || coords[1]>=dims[1]))
       ||(!periods[2] && (coords[2]<0 || coords[2]>=dims[2])))
      {
        // Unfortunately, the MPI standard for MPI_Cart_rank
        // was changed so that specifying coordinates outside
        // non-periodic dimensions results in an error;
        // therefore we now must take special care to set
        // these neighbors to null. Grrr.
        rank = MPI_PROC_NULL;
      }
      else
      {
        MPI_Cart_rank(CART_COMM, coords, &rank);
      }
      set_neighbor_rank(i0,i1,i2, rank);
    }
    // check that face neighbors are correct
    assert_eq(getLeftNeighbor(0),get_neighbor_rank(-1,0,0));
    assert_eq(getLeftNeighbor(1),get_neighbor_rank(0,-1,0));
    assert_eq(getLeftNeighbor(2),get_neighbor_rank(0,0,-1));
    assert_eq(getRghtNeighbor(0),get_neighbor_rank(1,0,0));
    assert_eq(getRghtNeighbor(1),get_neighbor_rank(0,1,0));
    assert_eq(getRghtNeighbor(2),get_neighbor_rank(0,0,1));
  }
  else {
    // previous check that nprocs = XLEN*YLEN*ZLEN should prevent reaching this line.
    eprintf("A process is thrown away from the new topology for fields.");
  }
  // This seems to be assumed elsewhere in the code.
  // We need to eliminate this assumption.
  // This becomes important if we want
  // a running program to have more than one
  // mesh and processor topology.
  // The MPI rank is for system-level code, e.g.
  // to identify the process from which debug is coming.
  // The cartesian rank is for application-level code.
  assert_eq(cartesian_rank, MPIdata::get_rank());

  for(int i=0;i<3;i++)
  {
    _isPeriodicLower[i] = periods[i] && isLower(i);
    _isPeriodicUpper[i] = periods[i] && isUpper(i);
  }

  for(int i=0;i<3;i++)
  {
    _noLeftNeighbor[i] = getLeftNeighbor(i) == MPI_PROC_NULL;
    _noRghtNeighbor[i] = getRghtNeighbor(i) == MPI_PROC_NULL;
  }

  _isBoundaryProcess = 
    noLeftNeighbor(0) ||
    noLeftNeighbor(1) ||
    noLeftNeighbor(2) ||
    noRghtNeighbor(0) ||
    noRghtNeighbor(1) ||
    noRghtNeighbor(2);
}
/** destructor */
VCtopology3D::~VCtopology3D() {

}
/** print topology info */
void VCtopology3D::Print() {
  cout << endl;
  cout << "Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Processors grid: "
    << dims[0] << "x" << dims[1] << "x" << dims[2] << endl;
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
  cout << "Neighbors: " << endl
       << "  xLeft = " << getLeftNeighbor(0) << endl
       << "  yLeft = " << getLeftNeighbor(1) << endl
       << "  zLeft = " << getLeftNeighbor(2) << endl
       << "  xRight = " << getRghtNeighbor(0) << endl
       << "  yRight = " << getRghtNeighbor(1) << endl
       << "  zRight = " << getRghtNeighbor(2) << endl;
  cout << endl;
}
