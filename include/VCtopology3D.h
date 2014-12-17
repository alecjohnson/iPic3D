/***************************************************************************
  VCtopology3D.h  -  a 3D Virtual cartesian topology
  A virtual topology is a mechanism for naming the processes
  in a communicator in a way that fits the communication
  pattern better. Since our processes will communicate mainly
  with the nearest neighbours after the fashion of a two-dimensional
  grid, we create a virtual topology to reflect this fact
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL Luveun
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef VCtopology3D_H
#define VCtopology3D_H

#include "mpi.h"

/**
 *  
 * Virtual cartesian topology
 * A virtual topology is a mechanism for naming the processes
 * in a communicator in a way that fits the communication
 * pattern better. Since our processes will communicate mainly
 * with the nearest neighbours after the fashion of a two-dimensional
 * grid, we create a virtual topology to reflect this fact
 * @version 2.0
 */
// 27 neighbors include self
const int NUM_COMM_NEIGHBORS=27;

class Collective;

class VCtopology3D //:public VirtualTopology3D
{
public:
  /** constructor: Define topology parameters: dimension, domain decomposition,... */
  VCtopology3D(const Collective& col);
  /** destructor */
  ~VCtopology3D();
  /** Find the neighbors in the new communicator  */
  void setup_vctopology(MPI_Comm comm_old);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();

  int get_dims(int i)const{return dims[i];}
  int getXLEN()const{ return dims[0]; }
  int getYLEN()const{ return dims[1]; }
  int getZLEN()const{ return dims[2]; }
  int getNprocs()const{ return (nprocs); }
  bool getPERIODICX()const{ return (periods[0]); }
  bool getPERIODICY()const{ return (periods[1]); }
  bool getPERIODICZ()const{ return (periods[2]); }

  // legacy names
  //
  int getCartesian_rank()const{ return (cartesian_rank); }
  int getXleft_neighbor() const{ return getLeftNeighbor(0); }
  int getYleft_neighbor() const{ return getLeftNeighbor(1); }
  int getZleft_neighbor() const{ return getLeftNeighbor(2); }
  int getXright_neighbor()const{ return getRghtNeighbor(0); }
  int getYright_neighbor()const{ return getRghtNeighbor(1); }
  int getZright_neighbor()const{ return getRghtNeighbor(2); }

  // new interface
  //
  bool is_rank0()const{ return (0==cartesian_rank); }
  int get_rank()const{ return (cartesian_rank); }
  int getLeftNeighbor(int i)const{return _leftNeighbor[i];}
  int getRghtNeighbor(int i)const{return _rghtNeighbor[i];}
  int getXleft()const{ return getLeftNeighbor(0); }
  int getYleft()const{ return getLeftNeighbor(1); }
  int getZleft()const{ return getLeftNeighbor(2); }
  int getXrght()const{ return getRghtNeighbor(0); }
  int getYrght()const{ return getRghtNeighbor(1); }
  int getZrght()const{ return getRghtNeighbor(2); }

  int get_neighbor_rank(int i, int j, int k)const
  {
    return neighbors3D[i+1][j+1][k+1];
  }
  int get_neighbor_rank(const int dirs[3])const
  {
    return neighbors3D[dirs[0]+1][dirs[1]+1][dirs[2]+1];
  }
  // forward map from {-1,0,1}x{-1,0,1}x{-1,0,1} to Z_{27}
  // based on row-major order:
  //
  //   00 01 02   09 10 11   18 19 20
  //   03 04 05   12 13 14   21 22 23
  //   06 07 08   15 16 17   24 25 26
  //
  static int get_index_for_direction(const int dirs[3])
  { return 3*(3*dirs[0]+dirs[1])+dirs[2]+13;}
  // inverse map from Z_{27} to {-1,0,1}x{-1,0,1}x{-1,0,1}
  static void set_direction_for_index(int dirs[3], int i)
  {
    dirs[2] = i%3-1; i/=3;
    dirs[1] = i%3-1; i/=3;
    dirs[0] = i%3-1;
  }
  static int get_negated_direction_tag(int di)
  {
    int dirs[3];
    set_direction_for_index(dirs,di);
    // flip the direction coordinates
    for(int i=0;i<3;i++)
    {
      dirs[i] *= -1;
    }
    return get_index_for_direction(dirs);
  }
  // flip the direction coordinate
  // in any direction that is pinned
  int get_pinned_direction_tag(int di)const
  {
    int dirs[3];
    set_direction_for_index(dirs,di);

    for(int i=0;i<3;i++)
    {
      const int coords = coordinates[i]+dirs[i];
      // pinning non-periodic coordinates
      if(!periods[i])
      {
        if(coords < 0 || coords >= dims[i])
        {
          dirs[i] *= -1;
        }
      }
    }
    return get_index_for_direction(dirs);
  }
  // (could precompute this on startup)
  // returns rank of neighbor in direction di
  int get_valid_neighbor_rank(int di)const
  {
    int dirs[3];
    set_direction_for_index(dirs,di);
    // "pin" displacement so that it points to the correct
    // distination after pinning
    for(int i=0;i<3;i++)
    {
      const int coords = coordinates[i]+dirs[i];
      // pin non-periodic coordinates
      if(!periods[i])
      {
        if((coords < 0) || (coords >= dims[i]))
        {
          dirs[i] = 0;
        }
      }
    }
    // return neighbor rank for pinned displacement
    return get_neighbor_rank(dirs);
  }

  bool isPeriodicLower(int i)const{ return _isPeriodicLower[i]; }
  bool isPeriodicUpper(int i)const{ return _isPeriodicUpper[i]; }
  //bool isPeriodicXlower()const{ return _isPeriodicLower[0]; }
  //bool isPeriodicXupper()const{ return _isPeriodicUpper[0]; }
  //bool isPeriodicYlower()const{ return _isPeriodicLower[1]; }
  //bool isPeriodicYupper()const{ return _isPeriodicUpper[1]; }
  //bool isPeriodicZlower()const{ return _isPeriodicLower[2]; }
  //bool isPeriodicZupper()const{ return _isPeriodicUpper[2]; }

  bool noLeftNeighbor(int i)const{return _noLeftNeighbor[i];}
  bool noRghtNeighbor(int i)const{return _noRghtNeighbor[i];}
  bool noXleftNeighbor()const{ return _noLeftNeighbor[0]; }
  bool noYleftNeighbor()const{ return _noLeftNeighbor[1]; }
  bool noZleftNeighbor()const{ return _noLeftNeighbor[2]; }
  bool noXrghtNeighbor()const{ return _noRghtNeighbor[0]; }
  bool noYrghtNeighbor()const{ return _noRghtNeighbor[1]; }
  bool noZrghtNeighbor()const{ return _noRghtNeighbor[2]; }

  bool hasXleftNeighbor()const{ return !_noLeftNeighbor[0]; }
  bool hasYleftNeighbor()const{ return !_noLeftNeighbor[1]; }
  bool hasZleftNeighbor()const{ return !_noLeftNeighbor[2]; }
  bool hasXrghtNeighbor()const{ return !_noRghtNeighbor[0]; }
  bool hasYrghtNeighbor()const{ return !_noRghtNeighbor[1]; }
  bool hasZrghtNeighbor()const{ return !_noRghtNeighbor[2]; }

  bool isBoundaryProcess()const{ return _isBoundaryProcess; }

  bool isLower(int i)const{ return coordinates[i]==0; }
  bool isUpper(int i)const{ return coordinates[i]==dims[i]-1; }

  bool getcVERBOSE()const{ return (cVERBOSE); }
  int getCoordinates(int dir)const{ return (coordinates[dir]); }
  const int *getCoordinates()const{ return (coordinates); }
  const int *getDims()const{ return dims; }
  //const int *getDivisions()const{ return getDims(); } // old name
  int getPeriods(int dir)const{ return (periods[dir]); }
  MPI_Comm getComm()const{ return (CART_COMM); }

private:

  void set_neighbor_rank(int i, int j, int k, int rank)
  {
    neighbors3D[i+1][j+1][k+1]=rank;
  }

private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;
  /** New communicator with virtual cartesian topology for Particles*/
  //MPI_Comm CART_COMM_P;
  /** MPI status during sending and receiving communication */
  MPI_Status status;
  /** Direction X for shift MPI_Cart_Shift*/
  int XDIR;
  /** Direction Y for shift MPI_Cart_Shift*/
  int YDIR;
  /** Direction Z for shift MPI_Cart_Shift*/
  int ZDIR;
  /** RIGHT = +    upwards   shift */
  int RIGHT;
  /** LEFT  = -    downwards shift */
  int LEFT;
  /** dimension of virtual topology */
  int PROCDIM;
  /** nprocs = number of processors */
  int nprocs;
  /** rank may be reordered     */
  int reorder;
  /** arrays for Create_Cart_create  */
  int dims[3]; // i.e. divisions
  /** periodicity */
  int periods[3];
  /** coordinates on processors grid */
  int coordinates[3];
  /** cartesian rank */
  int cartesian_rank;
  //
  // ranks of neighboring processes
  //
  int _leftNeighbor[3];
  int _rghtNeighbor[3];
  int neighbors3D[3][3][3];
  
  /** indicators of whether this is a periodic boundary */
  bool _isPeriodicLower[3];
  bool _isPeriodicUpper[3];

  /** indicators of whether this lacks a neighbor */
  bool _noLeftNeighbor[3];
  bool _noRghtNeighbor[3];

  int _isBoundaryProcess;

  /** if cVERBOSE == true, print to the screen all the comunication */
  bool cVERBOSE;
};

typedef VCtopology3D VirtualTopology3D;

#endif
