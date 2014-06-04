/*******************************************************************************************
  Particles3Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta.
 ********************************************************************************************/

#include <mpi.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include "asserts.h"
#include <algorithm> // for swap, std::max
#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "CollectiveIO.h"
#include "Collective.h"
#include "ComParticles3D.h"
#include "Alloc.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "Field.h"
#include "MPIdata.h"
#include "ompdefs.h"
#include "ipicmath.h"
#include "ipicdefs.h"

#include "Particle.h"
#include "Particles3Dcomm.h"
#include "Parameters.h"

#include "hdf5.h"
//#include <vector>
//#include <complex>
#include "debug.h"
#include "TimeTasks.h"

using std::cout;
using std::endl;

/**
 * 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** deallocate particles */
Particles3Dcomm::~Particles3Dcomm() {
  // extra xavg for sort
  MPI_Comm_free(&mpi_comm);
  delete numpcls_in_bucket;
  delete numpcls_in_bucket_now;
  delete bucket_offset;
}
/** constructor for a single species*/
// was Particles3Dcomm::allocate()
Particles3Dcomm::Particles3Dcomm(
  int species_number,
  CollectiveIO * col,
  VirtualTopology3D * vct_,
  Grid * grid_)
 :
  ns(species_number),
  vct(vct_),
  grid(grid_),
  pclIDgenerator(),
  particleType(ParticleType::AoS)
{
  // communicators for particles
  //
  MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
  // X direction
  sendXleft.init(Connection::null2self(vct->getXleft(),Connection::XDN,mpi_comm));
  sendXrght.init(Connection::null2self(vct->getXrght(),Connection::XUP,mpi_comm));
  recvXleft.init(Connection::null2self(vct->getXleft(),Connection::XUP,mpi_comm));
  recvXrght.init(Connection::null2self(vct->getXrght(),Connection::XDN,mpi_comm));
  // Y                                                                         
  sendYleft.init(Connection::null2self(vct->getYleft(),Connection::YDN,mpi_comm));
  sendYrght.init(Connection::null2self(vct->getYrght(),Connection::YUP,mpi_comm));
  recvYleft.init(Connection::null2self(vct->getYleft(),Connection::YUP,mpi_comm));
  recvYrght.init(Connection::null2self(vct->getYrght(),Connection::YDN,mpi_comm));
  // Z                                                                         
  sendZleft.init(Connection::null2self(vct->getZleft(),Connection::ZDN,mpi_comm));
  sendZrght.init(Connection::null2self(vct->getZrght(),Connection::ZUP,mpi_comm));
  recvZleft.init(Connection::null2self(vct->getZleft(),Connection::ZUP,mpi_comm));
  recvZrght.init(Connection::null2self(vct->getZrght(),Connection::ZDN,mpi_comm));

  recvXleft.post_recvs();
  recvXrght.post_recvs();
  recvYleft.post_recvs();
  recvYrght.post_recvs();
  recvZleft.post_recvs();
  recvZrght.post_recvs();

  // info from collectiveIO
  //
  npcel = col->getNpcel(get_species_num());
  npcelx = col->getNpcelx(get_species_num());
  npcely = col->getNpcely(get_species_num());
  npcelz = col->getNpcelz(get_species_num());
  //
  qom = col->getQOM(get_species_num());
  uth = col->getUth(get_species_num());
  vth = col->getVth(get_species_num());
  wth = col->getWth(get_species_num());
  u0 = col->getU0(get_species_num());
  v0 = col->getV0(get_species_num());
  w0 = col->getW0(get_species_num());
  dt = col->getDt();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  delta = col->getDelta();
  TrackParticleID = col->getTrackParticleID(get_species_num());
  c = col->getC();
  // info for mover
  NiterMover = col->getNiterMover();
  // velocity of the injection from the wall
  Vinj = col->getVinj();
  Ninj = col->getRHOinject(get_species_num());
  //
  // boundary condition for particles
  //
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceZright = col->getBcPfaceZright();
  bcPfaceZleft = col->getBcPfaceZleft();

  // info from Grid
  //
  xstart = grid->getXstart();
  xend = grid->getXend();
  ystart = grid->getYstart();
  yend = grid->getYend();
  zstart = grid->getZstart();
  zend = grid->getZend();
  //
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  inv_dx = 1/dx;
  inv_dy = 1/dy;
  inv_dz = 1/dz;
  //
  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();
  nxc = grid->getNXC();
  nyc = grid->getNYC();
  nzc = grid->getNZC();
  assert_eq(nxc,nxn-1);
  assert_eq(nyc,nyn-1);
  assert_eq(nzc,nzn-1);
  invVOL = grid->getInvVOL();

  // info from VirtualTopology3D
  //
  cVERBOSE = vct->getcVERBOSE();

  /////////////////////////////////
  // preallocate space in arrays //
  /////////////////////////////////
  //
  // determine number of particles to preallocate for this process.
  //
  // determine number of cells in this process
  //
  // we calculate in double precision to guard against overflow
  double dNp = double(grid->get_num_cells_rr())*col->getNpcel(species_number);
  double dNpmax = dNp * col->getNpMaxNpRatio();
  // ensure that particle index will not overflow 32-bit
  // representation as long as dmaxnop is respected.
  assert_le(dNpmax,double(INT_MAX));
  const int nop = dNp;
  // initialize particle ID generator based on number of particles
  // that will initially be produced.
  pclIDgenerator.reserve_num_particles(nop);
  // initialize each process with capacity for some extra particles
  const int initial_capacity = roundup_to_multiple(nop*1.2,DVECWIDTH);
  //
  // SoA particle representation
  //
  // velocities
  u.reserve(initial_capacity);
  v.reserve(initial_capacity);
  w.reserve(initial_capacity);
  // charge
  q.reserve(initial_capacity);
  // positions
  x.reserve(initial_capacity);
  y.reserve(initial_capacity);
  z.reserve(initial_capacity);
  // subcycle time
  t.reserve(initial_capacity);
  //
  // AoS particle representation
  //
  _pcls.reserve(initial_capacity);
  particleType = ParticleType::AoS; // canonical representation

  //
  // allocate arrays for sorting particles
  //
  numpcls_in_bucket = new array3_int(nxc,nyc,nzc);
  numpcls_in_bucket_now = new array3_int(nxc,nyc,nzc);
  bucket_offset = new array3_int(nxc,nyc,nzc);
  
  assert_eq(sizeof(SpeciesParticle),64);

  // if RESTART is true initialize the particle in allocate method
  restart = col->getRestart_status();
  if (restart != 0) {
    if (vct->getCartesian_rank() == 0 && get_species_num() == 0)
      cout << "LOADING PARTICLES FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[1];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << name_file << endl;
      cout << "RESTART NOT POSSIBLE" << endl;
    }

    stringstream species_name;
    species_name << get_species_num();
    // the cycle of the last restart is set to 0
    string name_dataset = "/particles/species_" + species_name.str() + "/x/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id); /* dataspace handle */
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    // get how many particles there are on this processor for this species
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    const int nop = dims_out[0]; // number of particles in this process
    // prepare arrays to receive particles
    particleType = ParticleType::SoA;
    resize_SoA(nop);
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x[0]);
    // close the data set
    status = H5Dclose(dataset_id);

    // get y
    name_dataset = "/particles/species_" + species_name.str() + "/y/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y[0]);
    status = H5Dclose(dataset_id);

    // get z
    name_dataset = "/particles/species_" + species_name.str() + "/z/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z[0]);
    status = H5Dclose(dataset_id);

    // get u
    name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    status = H5Dclose(dataset_id);
    // get v
    name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    status = H5Dclose(dataset_id);
    // get w
    name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &w[0]);
    status = H5Dclose(dataset_id);
    // get q
    name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &q[0]);
    status = H5Dclose(dataset_id);
    // ID 
    //if (TrackParticleID) {
    //  // herr_t (*old_func)(void*); // HDF 1.6
    //  H5E_auto2_t old_func;      // HDF 1.8.8
    //  void *old_client_data;
    //  H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8.8
    //  /* Turn off error handling */
    //  // H5Eset_auto(NULL, NULL); // HDF 1.6
    //  H5Eset_auto2(H5E_DEFAULT, 0, 0); // HDF 1.8
    //  name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
    //  dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    //
    //  // H5Eset_auto(old_func, old_client_data); // HDF 1.6
    //  H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    //  if (dataset_id > 0)
    //    status = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParticleID);
    //  //else {
    //  //  for (int counter = 0; counter < nop; counter++)
    //  //    fetch_ParticleID(counter) = particleIDgenerator.get_ID();
    //  //}
    //}
    // close the hdf file
    status = H5Fclose(file_id);
    convertParticlesToAoS();
  }
}

// pad capacities so that aligned vectorization
// does not result in an array overrun.
//
// this should usually be cheap (a no-op)
//
void Particles3Dcomm::pad_capacities()
{
  _pcls.reserve(roundup_to_multiple(_pcls.size(),DVECWIDTH));
  u.reserve(roundup_to_multiple(u.size(),DVECWIDTH));
  v.reserve(roundup_to_multiple(v.size(),DVECWIDTH));
  w.reserve(roundup_to_multiple(w.size(),DVECWIDTH));
  q.reserve(roundup_to_multiple(q.size(),DVECWIDTH));
  x.reserve(roundup_to_multiple(x.size(),DVECWIDTH));
  y.reserve(roundup_to_multiple(y.size(),DVECWIDTH));
  z.reserve(roundup_to_multiple(z.size(),DVECWIDTH));
  t.reserve(roundup_to_multiple(t.size(),DVECWIDTH));
}

void Particles3Dcomm::resize_AoS(int nop)
{
  const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
  _pcls.reserve(padded_nop);
  _pcls.resize(nop);
}

void Particles3Dcomm::resize_SoA(int nop)
{
  //
  // allocate space for particles including padding
  //
  const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
  u.reserve(padded_nop);
  v.reserve(padded_nop);
  w.reserve(padded_nop);
  q.reserve(padded_nop);
  x.reserve(padded_nop);
  y.reserve(padded_nop);
  z.reserve(padded_nop);
  t.reserve(padded_nop);
  //
  // define size of particle data
  //
  u.resize(nop);
  v.resize(nop);
  w.resize(nop);
  q.resize(nop);
  x.resize(nop);
  y.resize(nop);
  z.resize(nop);
  t.resize(nop);
}
// A much faster version of this is at EMfields3D::sumMoments
//
//void Particles3Dcomm::interpP2G(Field * EMf)
//{
//  const double inv_dx = 1.0 / dx;
//  const double inv_dy = 1.0 / dy;
//  const double inv_dz = 1.0 / dz;
//  const double nxn = grid->getNXN();
//  const double nyn = grid->getNYN();
//  const double nzn = grid->getNZN();
//  // assert_le(nop,(long long)INT_MAX); // else would need to use long long
//  // to make memory use scale to a large number of threads we
//  // could first apply an efficient parallel sorting algorithm
//  // to the particles and then accumulate moments in smaller
//  // subarrays.
//  {
//    for (int i = 0; i < nop; i++)
//    {
//      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
//      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
//      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
//      double temp[2][2][2];
//      double xi[2], eta[2], zeta[2];
//      xi[0] = x[i] - grid->getXN(ix - 1, iy, iz);
//      eta[0] = y[i] - grid->getYN(ix, iy - 1, iz);
//      zeta[0] = z[i] - grid->getZN(ix, iy, iz - 1);
//      xi[1] = grid->getXN(ix, iy, iz) - x[i];
//      eta[1] = grid->getYN(ix, iy, iz) - y[i];
//      zeta[1] = grid->getZN(ix, iy, iz) - z[i];
//      double weight[2][2][2];
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++) {
//            weight[ii][jj][kk] = q[i] * xi[ii] * eta[jj] * zeta[kk] * invVOL;
//          }
//      // add charge density
//      EMf->addRho(weight, ix, iy, iz, ns);
//      // add current density - X
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = u[i] * weight[ii][jj][kk];
//      EMf->addJx(temp, ix, iy, iz, ns);
//      // add current density - Y
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = v[i] * weight[ii][jj][kk];
//      EMf->addJy(temp, ix, iy, iz, ns);
//      // add current density - Z
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = w[i] * weight[ii][jj][kk];
//      EMf->addJz(temp, ix, iy, iz, ns);
//      // Pxx - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = u[i] * u[i] * weight[ii][jj][kk];
//      EMf->addPxx(temp, ix, iy, iz, ns);
//      // Pxy - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = u[i] * v[i] * weight[ii][jj][kk];
//      EMf->addPxy(temp, ix, iy, iz, ns);
//      // Pxz - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = u[i] * w[i] * weight[ii][jj][kk];
//      EMf->addPxz(temp, ix, iy, iz, ns);
//      // Pyy - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = v[i] * v[i] * weight[ii][jj][kk];
//      EMf->addPyy(temp, ix, iy, iz, ns);
//      // Pyz - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = v[i] * w[i] * weight[ii][jj][kk];
//      EMf->addPyz(temp, ix, iy, iz, ns);
//      // Pzz - add pressure tensor
//      for (int ii = 0; ii < 2; ii++)
//        for (int jj = 0; jj < 2; jj++)
//          for (int kk = 0; kk < 2; kk++)
//            temp[ii][jj][kk] = w[i] * w[i] * weight[ii][jj][kk];
//      EMf->addPzz(temp, ix, iy, iz, ns);
//    }
//  }
//  // communicate contribution from ghost cells 
//  EMf->communicateGhostP2G(ns, 0, 0, 0, 0, vct);
//}

// boundary conditions that must be applied to each dimension
//
// It would be better to call an SoA version of this method in
// the mover where particles are still arranged in SoA format,
// though in that case we would not be able to restrict the
// application of boundary conditions to boundary processes
// unless we also restrict the distance moved by a particle in
// each (subcycle) step.
//
//inline void Particles3Dcomm::apply_boundary_conditions(
//  SpeciesParticle& pcl,
//  bool isBoundaryProcess,
//  bool noXlowerNeighbor, bool noXupperNeighbor,
//  bool noYlowerNeighbor, bool noYupperNeighbor,
//  bool noZlowerNeighbor, bool noZupperNeighbor)
//{
//  if(isBoundaryProcess)
//  {
//    double& x = pcl.fetch_x();
//    double& y = pcl.fetch_y();
//    double& z = pcl.fetch_z();
//    double& u = pcl.fetch_u();
//    double& v = pcl.fetch_v();
//    double& w = pcl.fetch_w();
//    if (noXlowerNeighbor && pcl.get_x() < 0)
//      BCpclLeft(x,u,v,w,Lx,uth,vth,wth,bcPfaceXleft);
//    else if (noXupperNeighbor && pcl.get_x() > Lx)
//      BCpclRight(x,u,v,w,Lx,uth,vth,wth,bcPfaceXright);
//    if (noYlowerNeighbor && pcl.get_y() < 0)
//      BCpclLeft(y,u,v,w,Ly,uth,vth,wth,bcPfaceYleft);
//    else if (noYupperNeighbor && pcl.get_y() > Ly)
//      BCpclRight(y,u,v,w,Ly,uth,vth,wth,bcPfaceYright);
//    if (noZlowerNeighbor && pcl.get_z() < 0)
//      BCpclLeft(z,u,v,w,Lz,uth,vth,wth,bcPfaceZleft);
//    else if (noZupperNeighbor && pcl.get_z() > Lz)
//      BCpclRight(z,u,v,w,Lz,uth,vth,wth,bcPfaceZright);
//  }
//}

// returns true if particle was sent
//
// should vectorize this by comparing position vectors
//
inline bool Particles3Dcomm::send_pcl_to_appropriate_buffer(
  SpeciesParticle& pcl, int count[6])
{
  int was_sent = true;
  // put particle in appropriate communication buffer if exiting
  if(pcl.get_x() < xstart)
  {
    sendXleft.send(pcl);
    count[0]++;
  }
  else if(pcl.get_x() > xend)
  {
    sendXrght.send(pcl);
    count[1]++;
  }
  else if(pcl.get_y() < ystart)
  {
    sendYleft.send(pcl);
    count[2]++;
  }
  else if(pcl.get_y() > yend)
  {
    sendYrght.send(pcl);
    count[3]++;
  }
  else if(pcl.get_z() < zstart)
  {
    sendZleft.send(pcl);
    count[4]++;
  }
  else if(pcl.get_z() > zend)
  {
    sendZrght.send(pcl);
    count[5]++;
  }
  else was_sent = false;

  return was_sent;
}

// flush sending particles.
//
void Particles3Dcomm::flush_send()
{
  sendXleft.send_complete();
  sendXrght.send_complete();
  sendYleft.send_complete();
  sendYrght.send_complete();
  sendZleft.send_complete();
  sendZrght.send_complete();
}

// receive, sort, and, as appropriate, resend incoming particles
//
// assumes that flush_send() has been called
//
int Particles3Dcomm::handle_received_particles()
{
  // we expect to receive at least one communication from every
  // communicator, so make sure that all receive buffers are clear
  // and waiting
  //
  recvXleft.recv_start(); recvXrght.recv_start();
  recvYleft.recv_start(); recvYrght.recv_start();
  recvZleft.recv_start(); recvZrght.recv_start();

  // make sure that current block in each sender is ready for sending
  //
  sendXleft.send_start(); sendXrght.send_start();
  sendYleft.send_start(); sendYrght.send_start();
  sendZleft.send_start(); sendZrght.send_start();

  const int num_recv_buffers = 6;
  // determine the periodicity shift for each incoming buffer
  bool apply_shift[num_recv_buffers] =
  {
    vct->isPeriodicXlower(), vct->isPeriodicXupper(),
    vct->isPeriodicYlower(), vct->isPeriodicYupper(),
    vct->isPeriodicZlower(), vct->isPeriodicZupper()
  };
  bool apply_BCs[num_recv_buffers] =
  {
    vct->noXlowerNeighbor(), vct->noXupperNeighbor(),
    vct->noYlowerNeighbor(), vct->noYupperNeighbor(),
    vct->noZlowerNeighbor(), vct->noZupperNeighbor()
  };
  // The documentation in the input file says that boundary conditions
  // are simply ignored in the periodic case, so I omit this check.
  //for(int i=0;i<6;i++)assert(!(apply_shift[i]&&apply_BCs[i]));

  int recv_count[6]={0,0,0,0,0,0};
  int send_count[6]={0,0,0,0,0,0};
  int num_pcls_recved = 0;
  int num_pcls_resent = 0;
  // receive incoming particles, 
  // immediately resending any exiting particles
  //
  MPI_Request recv_requests[num_recv_buffers] = 
  {
    recvXleft.get_curr_request(), recvXrght.get_curr_request(),
    recvYleft.get_curr_request(), recvYrght.get_curr_request(),
    recvZleft.get_curr_request(), recvZrght.get_curr_request()
  };
  BlockCommunicator<SpeciesParticle>* recvBuffArr[num_recv_buffers] =
  {
    &recvXleft, &recvXrght,
    &recvYleft, &recvYrght,
    &recvZleft, &recvZrght
  };

  assert(!recvXleft.comm_finished());
  assert(!recvXrght.comm_finished());
  assert(!recvYleft.comm_finished());
  assert(!recvYrght.comm_finished());
  assert(!recvZleft.comm_finished());
  assert(!recvZrght.comm_finished());

  // while there are still incoming particles
  // put them in the appropriate buffer
  //
  while(!(
    recvXleft.comm_finished() && recvXrght.comm_finished() &&
    recvYleft.comm_finished() && recvYrght.comm_finished() &&
    recvZleft.comm_finished() && recvZrght.comm_finished()))
  {
    int recv_index;
    MPI_Status recv_status;
    dprintf("waiting for receive request");
    MPI_Waitany(num_recv_buffers, recv_requests, &recv_index, &recv_status);
    if(recv_index==MPI_UNDEFINED)
      eprintf("recv_requests contains no active handles");
    assert_ge(recv_index,0);
    assert_lt(recv_index,num_recv_buffers);
    //
    // grab the received block of particles and process it
    //
    BlockCommunicator<SpeciesParticle>* recvBuff = recvBuffArr[recv_index];
    Block<SpeciesParticle>& recv_block
      = recvBuff->fetch_received_block(recv_status);

    // if appropriate apply periodicity shift for this block
    //
    // I prefer to use modulo rather than a simple shift.
    // modulo guarantees an upper bound on the number of
    // communications needed before the particles are in
    // their correct domain.  Otherwise, the existence of
    // a single "OMG" particle can result in an unbounded
    // number of communication cycles.
    //
    if(apply_shift[recv_index])
    {
      aligned_vector(SpeciesParticle)& pcl_list = recv_block.fetch_block();
      switch(recv_index)
      {
        default:
          invalid_value_error(recv_index);
        case 0:
        case 1:
          double Lxinv = 1/Lx;
          for(int pidx=0;pidx<pcl_list.size();pidx++)
          {
            double& x = pcl_list[pidx].fetch_x();
            x = modulo(x, Lx, Lxinv);
            // if(recv_index==0) x -= Lx; else x += Lx;
          }
          break;
        case 2:
        case 3:
          double Lyinv = 1/Ly;
          for(int pidx=0;pidx<pcl_list.size();pidx++)
          {
            double& y = pcl_list[pidx].fetch_y();
            y = modulo(y, Ly, Lyinv);
            // if(recv_index==2) y -= Ly; else y += Ly;
          }
          break;
        case 4:
        case 5:
          double Lzinv = 1/Lz;
          for(int pidx=0;pidx<pcl_list.size();pidx++)
          {
            double& z = pcl_list[pidx].fetch_z();
            z = modulo(z, Lz, Lzinv);
            // if(recv_index==4) z -= Lz; else z += Lz;
          }
          break;
      }
    }
    // if appropriate apply boundary conditions to this block
    else if(apply_BCs[recv_index])
    {
      aligned_vector(SpeciesParticle)& pcl_list = recv_block.fetch_block();
      switch(recv_index)
      {
        default:
          invalid_value_error(recv_index);
        case 0: assert(vct->noXlowerNeighbor());
          apply_Xleft_BC(pcl_list);
          break;
        case 1: assert(vct->noXupperNeighbor());
          apply_Xrght_BC(pcl_list);
          break;
        case 2: assert(vct->noYlowerNeighbor());
          apply_Yleft_BC(pcl_list);
          break;
        case 3: assert(vct->noYupperNeighbor());
          apply_Yrght_BC(pcl_list);
          break;
        case 4: assert(vct->noZlowerNeighbor());
          apply_Zleft_BC(pcl_list);
          break;
        case 5: assert(vct->noZupperNeighbor());
          apply_Zrght_BC(pcl_list);
          break;
      }
    }
    dprintf("received %d particles from direction %s",
      recv_block.size(),
      recvBuff->get_connection().tag_name());

    recv_count[recv_index]+=recv_block.size();
    num_pcls_recved += recv_block.size();
    // process each particle in the received block.
    //
    for(int pidx=0;pidx<recv_block.size();pidx++)
    {
      SpeciesParticle& pcl = recv_block[pidx];
      dprintf("received particle %d", int(pcl.get_t()));
      bool was_sent = send_pcl_to_appropriate_buffer(pcl, send_count);

      if(__builtin_expect(was_sent,false))
      {
        num_pcls_resent++;
      }
      else
      {
        // the particle belongs here, so put it in the
        // appropriate place. for now all particles are in a
        // single list, so we append to the list.
        _pcls.push_back(pcl);
      }
    }
    // release the block and update the receive request
    recvBuff->release_received_block();
    recv_requests[recv_index] = recvBuff->get_curr_request();
  }

  dprintf("recv_count: %d+%d+%d+%d+%d+%d=%d",
    recv_count[0], recv_count[1], recv_count[2],
    recv_count[3], recv_count[4], recv_count[5], num_pcls_recved);
  dprintf("send_count: %d+%d+%d+%d+%d+%d=%d",
    send_count[0], send_count[1], send_count[2],
    send_count[3], send_count[4], send_count[5], num_pcls_resent);
  // return the number of particles that were resent

  return num_pcls_resent;
}

static long long mpi_global_sum(int in)
{
  long long total;
  long long long_in = (long long)in;
  dprintf("calling MPI_Allreduce(%d,&total,1, ...)", long_in);
  MPI_Allreduce(&long_in, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}

// these methods should be made virtual
// so that the user can override boundary conditions.
//
void Particles3Dcomm::apply_Xleft_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceXleft)
  {
    default:
      unsupported_value_error(bcPfaceXleft);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        pcls[p].fetch_x() *= -1;
        pcls[p].fetch_u() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      // in this case it might be faster to convert to and
      // from SoA format, if calls to rand() can vectorize.
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        pcl.fetch_x() *= -1;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        u[0] = fabs(u[0]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}
void Particles3Dcomm::apply_Yleft_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceYleft)
  {
    default:
      unsupported_value_error(bcPfaceYleft);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        pcls[p].fetch_y() *= -1;
        pcls[p].fetch_v() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      // in this case it might be faster to convert to and
      // from SoA format, if calls to rand() can vectorize.
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        pcl.fetch_y() *= -1;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        u[1] = fabs(u[1]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}
void Particles3Dcomm::apply_Zleft_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceZleft)
  {
    default:
      unsupported_value_error(bcPfaceZleft);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        pcls[p].fetch_z() *= -1;
        pcls[p].fetch_w() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      // in this case it might be faster to convert to and
      // from SoA format, if calls to rand() can vectorize.
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        pcl.fetch_z() *= -1;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        u[2] = fabs(u[2]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}
void Particles3Dcomm::apply_Xrght_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceXright)
  {
    default:
      unsupported_value_error(bcPfaceXright);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        double& x = pcls[p].fetch_x();
        x = 2*Lx - x;
        pcls[p].fetch_u() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      // in this case it might be faster to convert to and
      // from SoA format, if calls to rand() can vectorize.
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        double& x = pcl.fetch_x();
        x = 2*Lx - x;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        u[0] = -fabs(u[0]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}
void Particles3Dcomm::apply_Yrght_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceYright)
  {
    default:
      unsupported_value_error(bcPfaceYright);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        double& y = pcls[p].fetch_y();
        y = 2*Ly - y;
        pcls[p].fetch_v() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        double& y = pcl.fetch_y();
        y = 2*Ly - y;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        v[0] = -fabs(v[0]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}
void Particles3Dcomm::apply_Zrght_BC(
  aligned_vector(SpeciesParticle)& pcls)
{
  switch(bcPfaceZright)
  {
    default:
      unsupported_value_error(bcPfaceZright);
    case BCparticles::PERFECT_MIRROR:
      for(int p=0;p<pcls.size();p++)
      {
        double& z = pcls[p].fetch_z();
        z = 2*Lz - z;
        pcls[p].fetch_w() *= -1;
      }
      break;
    case BCparticles::REEMISSION:
      for(int p=0;p<pcls.size();p++)
      {
        SpeciesParticle& pcl = pcls[p];
        double& z = pcl.fetch_z();
        z = 2*Lz - z;
        double u[3];
        sample_maxwellian(u[0],u[1],u[2], uth,vth,wth);
        w[0] = -fabs(w[0]);
        pcl.set_u(u);
      }
      break;
    case BCparticles::EXIT:
      pcls.clear();
      break;
  }
}

// exchange particles with neighboring processors
//
// sent particles are deleted from _pcls.
// holes are filled with particles from end.
// then received particles are appended to end.
// returns number of unsent particles, i.e.,
// returns index of first received particle
//
int Particles3Dcomm::communicate_particles()
{
  const int num_ids = 1;
  longid id_list[num_ids] = {0};
  timeTasks_set_communicating(); // communicating until end of scope
  convertParticlesToAoS();
  print_pcls(_pcls,ns,id_list, num_ids);

  // activate receiving
  //
  recvXleft.recv_start(); recvXrght.recv_start();
  recvYleft.recv_start(); recvYrght.recv_start();
  recvZleft.recv_start(); recvZrght.recv_start();

  // make sure that current block in each sender is ready for sending
  //
  sendXleft.send_start(); sendXrght.send_start();
  sendYleft.send_start(); sendYrght.send_start();
  sendZleft.send_start(); sendZrght.send_start();

  int send_count[6]={0,0,0,0,0,0};
  const int orig_size = _pcls.size();
  int np_current = 0;
  while(np_current < _pcls.size())
  {
    SpeciesParticle& pcl = _pcls[np_current];
    // if the particle is exiting, put it in the appropriate send bin;
    // this could be done at conclusion of push after particles are
    // converted to AoS format in order to overlap communication
    // with computation.
    bool was_sent = send_pcl_to_appropriate_buffer(pcl,send_count);

    // fill in hole; for the sake of data pipelining could change
    // this to make a list of holes and then go back and fill
    // them in, but will builtin_expect also allow efficient
    // pipelining?  Or does the compiler generate instructions
    // that automatically adjust pipelining based on
    // accumulated statistical branching behavior?
    //
    // optimizer should assume that most particles are not sent
    if(__builtin_expect(was_sent,false))
    {
      //dprintf("sent particle %d", np_current);
      delete_particle(np_current);
    }
    else
    {
      np_current++;
    }
  }
  assert_eq(_pcls.size(),np_current);
  const int num_pcls_sent = orig_size - getNOP();
  //dprint(num_pcls_sent);
  dprintf("send_count = [%d,%d,%d,%d,%d,%d]",
    send_count[0], send_count[1], send_count[2],
    send_count[3], send_count[4], send_count[5]);
  dprintf("spec %d #pcls sent = %d", ns, num_pcls_sent);

  // most likely exactly three particle communications
  // will be needed, one for each dimension of space,
  // so we begin with three receives and thereafter
  // before each receive we do an all-reduce to check
  // if more particles actually need to be received.
  int num_pcls_resent;
  for(int i=0;i<3;i++)
  {
    // flush sending of particles
    flush_send();
    num_pcls_resent = handle_received_particles();
    //dprintf("spec %d #pcls resent = %d", ns, num_pcls_resent);
  }

  // continue receiving and resending incoming particles until
  // global all-reduce of num_pcls_resent is zero, indicating
  // that there are no more particles to be received.
  //
  long long total_num_pcls_resent = mpi_global_sum(num_pcls_resent);
  dprintf("spec %d pcls resent = %d, %d",
    ns, num_pcls_resent, total_num_pcls_resent);
  // the maximum number of neighbor communications needed to
  // put a particle in the correct mesh cell.
  int comm_limit = 3*std::max(vct->getXLEN(),
    std::max(vct->getYLEN(), vct->getYLEN()));
  int comm_idx=0;
  while(total_num_pcls_resent)
  {
    if(comm_idx>=comm_limit);
      eprintf("failed to complete particle communication"
        " within %d communications", comm_limit);

    // flush sending of particles
    flush_send();
    num_pcls_resent = handle_received_particles();

    total_num_pcls_resent = mpi_global_sum(num_pcls_resent);
    dprint(total_num_pcls_resent);
    comm_idx++;
  }
  //print_pcls(_pcls,ns,0,0);
}

/** return the Kinetic energy */
double Particles3Dcomm::getKe() {
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register int i = 0; i < _pcls.size(); i++)
  {
    SpeciesParticle& pcl = _pcls[i];
    const double u = pcl.get_u();
    const double v = pcl.get_v();
    const double w = pcl.get_w();
    const double q = pcl.get_q();
    localKe += .5*(q/qom)*(u*u + v*v + w*w);
  }
  MPI_Allreduce(&localKe, &totalKe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalKe);
}

/** return the total momentum */
//
// This is the sum over all particles of the magnitude of the
// momentum, which has no physical meaning that I can see.
// we should be summing each component of the momentum. -eaj
//
double Particles3Dcomm::getP() {
  double localP = 0.0;
  double totalP = 0.0;
  for (register int i = 0; i < _pcls.size(); i++)
  {
    SpeciesParticle& pcl = _pcls[i];
    const double u = pcl.get_u();
    const double v = pcl.get_v();
    const double w = pcl.get_w();
    const double q = pcl.get_q();
    localP += (q/qom)*sqrt(u*u + v*v + w*w);
  }
  MPI_Allreduce(&localP, &totalP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalP);
}

/** return the highest kinetic energy */
double Particles3Dcomm::getMaxVelocity() {
  double localVel = 0.0;
  double maxVel = 0.0;
  for (int i = 0; i < _pcls.size(); i++)
  {
    SpeciesParticle& pcl = _pcls[i];
    const double u = pcl.get_u();
    const double v = pcl.get_v();
    const double w = pcl.get_w();
    localVel = std::max(localVel, sqrt(u*u + v*v + w*w));
  }
  MPI_Allreduce(&localVel, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return (maxVel);
}


/** get energy spectrum */
//
// this ignores the weight of the charges. -eaj
//
long long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel) {
  long long *f = new long long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (int i = 0; i < _pcls.size(); i++) {
    SpeciesParticle& pcl = _pcls[i];
    const double u = pcl.get_u();
    const double v = pcl.get_v();
    const double w = pcl.get_w();
    Vel = sqrt(u*u + v*v + w*w);
    bin = int (floor(Vel / dv));
    if (bin >= nBins)
      f[nBins - 1] += 1;
    else
      f[bin] += 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, f, nBins, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  // This way of summing is very inefficient
  //{
  //  long long localN = 0;
  //  long long totalN = 0;
  //  for (int i = 0; i < nBins; i++) {
  //    localN = f[i];
  //    MPI_Allreduce(&localN, &totalN, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  //    f[i] = totalN;
  //  }
  //}
  return f;
}


/** print particles info */
void Particles3Dcomm::Print() const
{
  cout << endl;
  cout << "Number of Particles: " << _pcls.size() << endl;
  cout << "Subgrid (" << vct->getCoordinates(0) << "," << vct->getCoordinates(1) << "," << vct->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << get_species_num() << endl;
  for (int i = 0; i < _pcls.size(); i++)
  {
    const SpeciesParticle& pcl = _pcls[i];
    cout << "Particle #" << i << ":"
      << " x=" << pcl.get_x()
      << " y=" << pcl.get_y()
      << " z=" << pcl.get_z()
      << " u=" << pcl.get_u()
      << " v=" << pcl.get_v()
      << " w=" << pcl.get_w()
      << endl;
  }
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp()  const
{
  cout << endl;
  cout << "Number of Particles of species " << get_species_num() << ": " << getNOP() << endl;
  cout << "Subgrid (" << vct->getCoordinates(0) << "," << vct->getCoordinates(1) << "," << vct->getCoordinates(2) << ")" << endl;
  cout << endl;
}

/***** particle sorting routines *****/

void Particles3Dcomm::sort_particles_serial()
{
  switch(particleType)
  {
    case ParticleType::AoS:
      sort_particles_serial_AoS();
      break;
    case ParticleType::SoA:
      convertParticlesToAoS();
      sort_particles_serial_AoS();
      convertParticlesToSynched();
      break;
    default:
      unsupported_value_error(particleType);
  }
}

// need to sort and communicate particles after each iteration
void Particles3Dcomm::sort_particles_serial_AoS()
{
  convertParticlesToAoS();

  _pclstmp.reserve(_pcls.size());
  {
    numpcls_in_bucket->setall(0);
    // iterate through particles and count where they will go
    for (int pidx = 0; pidx < _pcls.size(); pidx++)
    {
      const SpeciesParticle& pcl = get_pcl(pidx);
      // get the cell indices of the particle
      int cx,cy,cz;
      grid->get_safe_cell_coordinates(cx,cy,cz,pcl.get_x(),pcl.get_y(),pcl.get_z());

      // increment the number of particles in bucket of this particle
      (*numpcls_in_bucket)[cx][cy][cz]++;
    }

    // compute prefix sum to determine initial position
    // of each bucket (could parallelize this)
    //
    int accpcls=0;
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    {
      (*bucket_offset)[cx][cy][cz] = accpcls;
      accpcls += (*numpcls_in_bucket)[cx][cy][cz];
    }
    assert_eq(accpcls,getNOP());

    numpcls_in_bucket_now->setall(0);
    // put the particles where they are supposed to go
    const int nop = getNOP();
    for (int pidx = 0; pidx < nop; pidx++)
    {
      const SpeciesParticle& pcl = get_pcl(pidx);
      // get the cell indices of the particle
      int cx,cy,cz;
      grid->get_safe_cell_coordinates(cx,cy,cz,pcl.get_x(),pcl.get_y(),pcl.get_z());

      // compute where the data should go
      const int numpcls_now = (*numpcls_in_bucket_now)[cx][cy][cz]++;
      const int outpidx = (*bucket_offset)[cx][cy][cz] + numpcls_now;
      assert_lt(outpidx, nop);
      assert_ge(outpidx, 0);
      assert_lt(pidx, nop);
      assert_ge(pidx, 0);

      // copy particle data to new location
      //
      _pclstmp[outpidx] = pcl;
    }
    // swap the tmp particle memory with the official particle memory
    {
      // if using accessors rather than transposition,
      // here I would need not only to swap the pointers but also
      // to swap all the accessors.
      //
      _pcls.swap(_pclstmp);
    }

    // check if the particles were sorted incorrectly
    if(true)
    {
      for(int cx=0;cx<nxc;cx++)
      for(int cy=0;cy<nyc;cy++)
      for(int cz=0;cz<nzc;cz++)
      {
        assert_eq((*numpcls_in_bucket_now)[cx][cy][cz], (*numpcls_in_bucket)[cx][cy][cz]);
      }
    }
  }
  // SoA particle representation is no longer valid
  particleType = ParticleType::AoS;
}

//void Particles3Dcomm::sort_particles_parallel(
//  double *xpos, double *ypos, double *zpos,
//  Grid * grid, VirtualTopology3D * vct)
//{
//  // should change this to first communicate particles so that
//  // they are in the correct process and all particles
//  // lie in this subdomain.
//
//  // count the number of particles to go in each bucket
//  numpcls_in_bucket.setall(0);
//  #pragma omp parallel
//  {
//    const int thread_num = omp_get_thread_num();
//    arr3_int numpcls_in_bucket_thr = fetch_numpcls_in_bucket_thr(thread_num);
//    numpcls_in_bucket_thr.setall(0);
//    // iterate through particles and count where they will go
//    #pragma omp for // nowait
//    for (int pidx = 0; pidx < nop; pidx++)
//    {
//      // get the cell indices of the particle
//      // (should change this to use xavg[pidx])
//      const pfloat xpos = xpos[pidx];
//      const pfloat ypos = ypos[pidx];
//      const pfloat zpos = zpos[pidx];
//      int cx,cy,cz;
//      get_safe_cell_for_pos(cx,cy,cz,xpos,ypos,zpos);
//
//      // need to allocate these
//      //
//      //xidx[pidx]=cx;
//      //yidx[pidx]=cy;
//      //zidx[pidx]=cz;
//
//      // increment the number of particles in bucket of this particle
//      numpcls_in_bucket_thr[cx][cy][cz]++;
//    }
//    // reduce the thread buckets into the main bucket
//    // #pragma omp critical (numpcls_in_bucket_reduction)
//    {
//      #pragma omp for collapse(2)
//      for(int cx=0;cx<nxc;cx++)
//      for(int cy=0;cy<nyc;cy++)
//      for(int th=0;th<num_threads;th++)
//      for(int cz=0;cz<nzc;cz++)
//      {
//        numpcls_in_bucket[cx][cy][cz]
//          += get_numpcls_in_bucket_thr(th)[cx][cy][cz];
//      }
//    }
//
//    // compute prefix sum to determine initial position
//    // of each bucket (could parallelize this)
//    //
//    int accpcls=0;
//    #pragma omp critical (bucket_offset_reduction)
//    for(int cx=0;cx<nxc;cx++)
//    for(int cy=0;cy<nyc;cy++)
//    for(int cz=0;cz<nzc;cz++)
//    {
//      bucket_offset[cx][cy][cz] = accpcls;
//      accpcls += numpcls_in_bucket[cx][cy][cz];
//    }
//
//    // cycle through the mesh cells mod 3
//    // (or mod(2*N+1), where N is number of mesh cells
//    // that a slow particle can move).
//    // This ensures that slow particles can be moved
//    // to their destinations without write conflicts
//    // among threads.  But what about cache contention?
//    //
//    for(int cxmod3=0; cxmod3<3; cxmod3++)
//    #pragma omp for collapse(2)
//    for(int cx=cxmod3; cx<nxc; cx+=3)
//    for(int cy=0; cy<nyc; cy++)
//    for(int cz=0; cz<nzc; cz++)
//    {
//      // put the slow particles where they are supposed to go and
//      // set aside the fast particles for separate processing.
//      // (to vectorize would need to sort separately in each
//      // dimension of space).
//      //
//      // problem: particles might have to be moved not because
//      // they are fast but because of an overall shift in the
//      // number of particles in a location, e.g. because of
//      // particles flowing in from a jet. Need a different
//      // approach, where memory is allocated for each cell.
//      _numpcls_in_bucket = numpcls_in_bucket[cx][cy][cz];
//      for(int pidx=bucket_offset[cx][cy][cz]; pidx<_numpcls_in_bucket; pidx++)
//      {
//        const int outcx = xidx[pidx];
//        const int outcy = yidx[pidx];
//        const int outcz = zidx[pidx];
//        const int cxlower = outcx <= 0 ? 0 : outcx-1;
//        const int cxupper = outcx >= (nxc-1) ? nxc-1 : outcx+1;
//        const int lowerindex = bucket_offset[cxlower][cylower][czlower];
//        const int upperoffset = bucket_offset[cxupper][cyupper][czupper];
//        const int upperindex = upperoffset + numpcls_in_bucket[outcx][outcy][outcz];
//        ...
//      }
//    }
//    // (1) put fast particles that must be moved more than one
//    // mesh cell at the end of the cell's list, and
//    // (2) put slow particles in the correct location
//
//    // count the number of particles that need to be moved
//    // more than one mesh cell and allocate a special buffer for them.
//    // (could change to count number of particles that need
//    // to move more than N mesh cells).
//    //
//    int numpcls_long_move_thr = 0;
//    #pragma omp for // nowait
//    for (int i = 0; i < nop; i++)
//    {
//      const int cx = xidx[pidx];
//      const int cy = yidx[pidx];
//      const int cz = zidx[pidx];
//
//      const int cxlower = cx <= 0 ? 0 : cx-1;
//      const int cxupper = cx >= (nxc-1) ? nxc-1 : cx+1;
//      const int lowerindex = bucket_offset[cxlower][cylower][czlower];
//      const int upperoffset = bucket_offset[cxupper][cyupper][czupper];
//      const int upperindex = upperoffset + numpcls_in_bucket[cx][cy][cz];
//      if(i < lowerindex || i > upperindex)
//      {
//        numpcls_long_move_thr++;
//      }
//    }
//  }
//}
//#endif

void Particles3Dcomm::copyParticlesToSoA()
{
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_SOA);
  dprintf("copying to struct of arrays");
  const int nop = _pcls.size();
  // create memory for SoA representation
  resize_SoA(nop);
 #ifndef __MIC__
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx++)
  {
    const SpeciesParticle& pcl = _pcls[pidx];
    u[pidx] = pcl.get_u();
    v[pidx] = pcl.get_v();
    w[pidx] = pcl.get_w();
    q[pidx] = pcl.get_q();
    x[pidx] = pcl.get_x();
    y[pidx] = pcl.get_y();
    z[pidx] = pcl.get_z();
    t[pidx] = pcl.get_t();
  }
 #else // __MIC__
  // rather than doing stride-8 scatter,
  // copy and transpose data 8 particles at a time
  assert_divides(8,u.capacity());
  for(int pidx=0; pidx<nop; pidx+=8)
  {
    (F64vec8*) SoAdata[8] = {&u[pidx],&v[pidx],&w[pidx],&q[pidx],
                             &x[pidx],&y[pidx],&z[pidx],&t[pidx]};
    F64vec8 (&AoSdata)[8] = *reinterpret_cast<double (*)F64vec8[8]>(&pcls[pidx]);
    transpose_8x8_double(AoSdata,SoAdata);
  }
 #endif // __MIC__
  particleType = ParticleType::synched;
}

void Particles3Dcomm::copyParticlesToAoS()
{
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_AOS);
  const int nop = u.size();
  resize_AoS(nop);
  dprintf("copying to array of structs");
 #ifndef __MIC__
  // use a simple stride-8 gather
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx++)
  {
    _pcls[pidx].set(
      u[pidx],v[pidx],w[pidx], q[pidx],
      x[pidx],y[pidx],z[pidx], t[pidx]);
  }
 #else // __MIC__
  // for efficiency, copy data 8 particles at a time,
  // transposing each block of particles
  assert_divides(8,_pcls.capacity());
  for(int pidx=0; pidx<nop; pidx+=8);
  {
    //double (&matrix)[8][8] = *(double (*)[8][8])(_pcls[pidx]);
    //double (&matrix)[8][8] = *reinterpret_cast<double (*)[8][8]>(&_pcls[pidx]);
    F64vec8 (AoSdata&) [8] = *reinterpret_cast<F64vec8(*)[8]>(&_pcls[pidx]);
    (F64vec8*) SoAdata[8] = {&u[pidx],&v[pidx],&w[pidx],&q[pidx],
                             &x[pidx],&y[pidx],&z[pidx],&t[pidx]};
    transpose_8x8_double(SoAdata, AoSdata);
  }
 #endif
  particleType = ParticleType::synched;
}

// synched AoS and SoA conceptually implies a write-lock
//
void Particles3Dcomm::convertParticlesToSynched()
{
  switch(particleType)
  {
    default:
      unsupported_value_error(particleType);
    case ParticleType::SoA:
      copyParticlesToAoS();
      break;
    case ParticleType::AoS:
      copyParticlesToSoA();
      break;
    case ParticleType::synched:
      break;
  }
  // this state conceptually implies a write-lock
  particleType = ParticleType::synched;
}

// defines AoS to be the authority
// (conceptually releasing any write-lock)
//
void Particles3Dcomm::convertParticlesToAoS()
{
  switch(particleType)
  {
    default:
      unsupported_value_error(particleType);
    case ParticleType::SoA:
      copyParticlesToAoS();
      break;
    case ParticleType::AoS:
    case ParticleType::synched:
      break;
  }
  particleType = ParticleType::AoS;
}

// check whether particles are SoA
bool Particles3Dcomm::particlesAreSoA()const
{
  switch(particleType)
  {
    default:
      unsupported_value_error(particleType);
    case ParticleType::AoS:
      return false;
      break;
    case ParticleType::SoA:
    case ParticleType::synched:
      break;
  }
  return true;
}

// defines SoA to be the authority
// (conceptually releasing any write-lock)
//
void Particles3Dcomm::convertParticlesToSoA()
{
  switch(particleType)
  {
    default:
      unsupported_value_error(particleType);
    case ParticleType::AoS:
      copyParticlesToSoA();
      break;
    case ParticleType::SoA:
    case ParticleType::synched:
      break;
  }
  particleType = ParticleType::SoA;
}

