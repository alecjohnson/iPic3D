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
  particleType(ParticleType::AoS)
{
  // info from collectiveIO
  //
  npcel = col->getNpcel(get_species_num());
  npcelx = col->getNpcelx(get_species_num());
  npcely = col->getNpcely(get_species_num());
  npcelz = col->getNpcelz(get_species_num());
  //
  // determine number of particles to preallocate for this process.
  //
  int nop = col->getNp(get_species_num()) / (vct->getNprocs());
  //int npmax = 2*nop;
  int npmax = col->getNpMax(get_species_num()) / (vct->getNprocs());
  //np_tot = col->getNp(get_species_num());
  npmax = roundup_to_multiple(npmax,AoS_PCLS_AT_A_TIME);
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
  // SoA particle representation
  //
  // velocities
  u.reserve(npmax);
  v.reserve(npmax);
  w.reserve(npmax);
  // charge
  q.reserve(npmax);
  // positions
  x.reserve(npmax);
  y.reserve(npmax);
  z.reserve(npmax);
  // subcycle time
  t.reserve(npmax);
  //
  // AoS particle representation
  //
  _pcls.reserve(npmax);

  // communicators for particles
  //
  const MPI_Comm comm = MPI_COMM_WORLD;
  // X direction
  sendXleft.init(Connection(vct->getXleft(),Connection::PARTICLE_DN,comm));
  sendXrght.init(Connection(vct->getXrght(),Connection::PARTICLE_UP,comm));
  recvXleft.init(Connection(vct->getXleft(),Connection::PARTICLE_UP,comm));
  recvXrght.init(Connection(vct->getXrght(),Connection::PARTICLE_DN,comm));
  // Y direction
  sendYleft.init(Connection(vct->getYleft(),Connection::PARTICLE_DN,comm));
  sendYrght.init(Connection(vct->getYrght(),Connection::PARTICLE_UP,comm));
  recvYleft.init(Connection(vct->getYleft(),Connection::PARTICLE_UP,comm));
  recvYrght.init(Connection(vct->getYrght(),Connection::PARTICLE_DN,comm));
  // Z direction
  sendZleft.init(Connection(vct->getZleft(),Connection::PARTICLE_DN,comm));
  sendZrght.init(Connection(vct->getZrght(),Connection::PARTICLE_UP,comm));
  recvZleft.init(Connection(vct->getZleft(),Connection::PARTICLE_UP,comm));
  recvZrght.init(Connection(vct->getZrght(),Connection::PARTICLE_DN,comm));
  //
  // allocate arrays for sorting particles
  //
  numpcls_in_bucket = new array3_int(nxc,nyc,nzc);
  numpcls_in_bucket_now = new array3_int(nxc,nyc,nzc);
  bucket_offset = new array3_int(nxc,nyc,nzc);
  
  if(Parameters::get_USING_AOS())
  {
    assert_eq(sizeof(SpeciesParticle),64);
  }

  //ParticleID = new Larray<longid>;

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
    const int nop = dims_out[0];          // this the number of particles on the processor!
    u.resize(nop);
    v.resize(nop);
    w.resize(nop);
    q.resize(nop);
    x.resize(nop);
    y.resize(nop);
    z.resize(nop);
    t.resize(nop);
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
  }

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
inline bool Particles3Dcomm::apply_boundary_conditions(
  SpeciesParticle& pcl,
  bool isBoundaryProcess,
  bool noXlowerNeighbor, bool noXupperNeighbor,
  bool noYlowerNeighbor, bool noYupperNeighbor,
  bool noZlowerNeighbor, bool noZupperNeighbor)
{
  if(isBoundaryProcess)
  {
    double& x = pcl.fetch_x();
    double& y = pcl.fetch_y();
    double& z = pcl.fetch_z();
    double& u = pcl.fetch_u();
    double& v = pcl.fetch_v();
    double& w = pcl.fetch_w();
    if (noXlowerNeighbor && pcl.get_x() < 0)
      BCpclLeft(x,u,v,w,Lx,uth,vth,wth,bcPfaceXleft);
    else if (noXupperNeighbor && pcl.get_x() > Lx)
      BCpclRight(x,u,v,w,Lx,uth,vth,wth,bcPfaceXright);
    if (noYlowerNeighbor && pcl.get_y() < 0)
      BCpclLeft(y,u,v,w,Ly,uth,vth,wth,bcPfaceYleft);
    else if (noYupperNeighbor && pcl.get_y() > Ly)
      BCpclRight(y,u,v,w,Ly,uth,vth,wth,bcPfaceYright);
    if (noZlowerNeighbor && pcl.get_z() < 0)
      BCpclLeft(z,u,v,w,Lz,uth,vth,wth,bcPfaceZleft);
    else if (noZupperNeighbor && pcl.get_z() > Lz)
      BCpclRight(z,u,v,w,Lz,uth,vth,wth,bcPfaceZright);
  }
}

// returns true if particle was sent
//
inline bool Particles3Dcomm::send_pcl_to_appropriate_buffer(
  SpeciesParticle& pcl,
  bool hasXlowerNeighbor, bool hasXupperNeighbor,
  bool hasYlowerNeighbor, bool hasYupperNeighbor,
  bool hasZlowerNeighbor, bool hasZupperNeighbor,
  bool isPeriodicXlower, bool isPeriodicXupper,
  bool isPeriodicYlower, bool isPeriodicYupper,
  bool isPeriodicZlower, bool isPeriodicZupper)
{
  bool was_sent = true;
  
  // put particle in appropriate communication buffer if exiting
  //
  // (should change to do this immediately after pushing it so that
  // only once pass through particles is necessary and so that
  // communication can overlap computation)
  //
  if(hasXlowerNeighbor && pcl.get_x() < xstart)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicXlower && pcl.get_x() < 0) pcl.fetch_x() += Lx;
    // put it in the communication buffer
    sendXleft.send(pcl);
  }
  else if(hasXupperNeighbor && pcl.get_x() > xend)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicXupper && pcl.get_x() > Lx) pcl.fetch_x() -= Lx;
    // put it in the communication buffer
    sendXleft.send(pcl);
  }
  else if(hasYlowerNeighbor && pcl.get_y() < ystart)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicYlower && pcl.get_y() < 0) pcl.fetch_y() += Ly;
    // put it in the communication buffer
    sendYleft.send(pcl);
  }
  else if(hasYupperNeighbor && pcl.get_y() > yend)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicYupper && pcl.get_y() > Ly) pcl.fetch_y() -= Ly;
    sendYrght.send(pcl);
  }
  else if(hasZlowerNeighbor && pcl.get_z() < zstart)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicZlower && pcl.get_z() < 0) pcl.fetch_z() += Lz;
    // put it in the communication buffer
    sendZleft.send(pcl);
  }
  else if(hasZupperNeighbor && pcl.get_z() > zend)
  {
    // handle periodic boundary conditions only when wrapping particles
    if(isPeriodicZupper && pcl.get_z() > Lz) pcl.fetch_z() -= Lz;
    // put it in the communication buffer
    sendZrght.send(pcl);
  }
  else {
    // particle is still in the domain, procede with the next particle
    was_sent = false;
  }
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
int Particles3Dcomm::handle_incoming_particles()
{
  int num_pcls_resent = 0;
  // receive incoming particles, 
  // immediately resending any exiting particles
  //
  const int num_recv_buffers = 6;
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

  const bool noXlowerNeighbor = vct->noXlowerNeighbor();
  const bool noXupperNeighbor = vct->noXupperNeighbor();
  const bool noYlowerNeighbor = vct->noYlowerNeighbor();
  const bool noYupperNeighbor = vct->noYupperNeighbor();
  const bool noZlowerNeighbor = vct->noZlowerNeighbor();
  const bool noZupperNeighbor = vct->noZupperNeighbor();
  const bool isBoundaryProcess = vct->isBoundaryProcess();
  const bool hasXlowerNeighbor = !vct->noXlowerNeighbor();
  const bool hasXupperNeighbor = !vct->noXupperNeighbor();
  const bool hasYlowerNeighbor = !vct->noYlowerNeighbor();
  const bool hasYupperNeighbor = !vct->noYupperNeighbor();
  const bool hasZlowerNeighbor = !vct->noZlowerNeighbor();
  const bool hasZupperNeighbor = !vct->noZupperNeighbor();
  const bool isPeriodicXlower = vct->isPeriodicXlower();
  const bool isPeriodicXupper = vct->isPeriodicXupper();
  const bool isPeriodicYlower = vct->isPeriodicYlower();
  const bool isPeriodicYupper = vct->isPeriodicYupper();
  const bool isPeriodicZlower = vct->isPeriodicZlower();
  const bool isPeriodicZupper = vct->isPeriodicZupper();

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
    MPI_Waitany(num_recv_buffers, recv_requests, &recv_index, &recv_status);
    if(recv_index==MPI_UNDEFINED)
      eprintf("recv_requests contains no active handles");
    assert_ge(recv_index,0);
    assert_lt(recv_index,num_recv_buffers);
    //
    // code specific to the receive buffer could be handled with a
    // switch(recv_index) code block.
    //
    // grab the received block of particles and process it
    //
    BlockCommunicator<SpeciesParticle>* recvBuff = recvBuffArr[recv_index];
    Block<Particle>& recv_block = recvBuff->fetch_received_block(recv_status);

    // process each particle in the received block.
    //
    for(int i=0;i<recv_block.size();i++)
    {
      SpeciesParticle& pcl = recv_block[i];
      apply_boundary_conditions(pcl,
        isBoundaryProcess,
        noXlowerNeighbor, noXupperNeighbor,
        noYlowerNeighbor, noYupperNeighbor,
        noZlowerNeighbor, noZupperNeighbor);
      bool was_sent = send_pcl_to_appropriate_buffer(pcl,
        isBoundaryProcess,
        hasXlowerNeighbor, hasXupperNeighbor,
        hasYlowerNeighbor, hasYupperNeighbor,
        hasZlowerNeighbor, hasZupperNeighbor,
        isPeriodicXlower, isPeriodicXupper,
        isPeriodicYlower, isPeriodicYupper,
        isPeriodicZlower, isPeriodicZupper)

      if(was_sent)
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
    recvBuff->release_received_block(recv_status);
    recv_requests[recv_index] = recvBuff->get_curr_request();
  }

  // return the number of particles that were resent
  return num_pcls_resent;
}

static long long mpi_global_sum(int in)
{
  long long total;
  long long long_in = long long(in);
  MPI_Allreduce(&long_in, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  RESOLVE: WHAT DOES // THE 1 MEAN?
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
  timeTasks_set_communicating(); // communicating until end of scope
  convertParticlesToAoS();

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

  const bool noXlowerNeighbor = ptVCT->noXlowerNeighbor();
  const bool noXupperNeighbor = ptVCT->noXupperNeighbor();
  const bool noYlowerNeighbor = ptVCT->noYlowerNeighbor();
  const bool noYupperNeighbor = ptVCT->noYupperNeighbor();
  const bool noZlowerNeighbor = ptVCT->noZlowerNeighbor();
  const bool noZupperNeighbor = ptVCT->noZupperNeighbor();
  const bool isBoundaryProcess = ptVCT->isBoundaryProcess();
  const bool hasXlowerNeighbor = !ptVCT->noXlowerNeighbor();
  const bool hasXupperNeighbor = !ptVCT->noXupperNeighbor();
  const bool hasYlowerNeighbor = !ptVCT->noYlowerNeighbor();
  const bool hasYupperNeighbor = !ptVCT->noYupperNeighbor();
  const bool hasZlowerNeighbor = !ptVCT->noZlowerNeighbor();
  const bool hasZupperNeighbor = !ptVCT->noZupperNeighbor();
  const bool isPeriodicXlower = ptVCT->isPeriodicXlower();
  const bool isPeriodicXupper = ptVCT->isPeriodicXupper();
  const bool isPeriodicYlower = ptVCT->isPeriodicYlower();
  const bool isPeriodicYupper = ptVCT->isPeriodicYupper();
  const bool isPeriodicZlower = ptVCT->isPeriodicZlower();
  const bool isPeriodicZupper = ptVCT->isPeriodicZupper();

  int np_current = 0;
  while(np_current < _pcls.size())
  {
    SpeciesParticle& pcl = _pcls[np_current];
    // should change to enforce boundary conditions at conclusion of push,
    // when particles are still in SoA format.
    apply_boundary_conditions(pcl,
      isBoundaryProcess,
      noXlowerNeighbor, noXupperNeighbor,
      noYlowerNeighbor, noYupperNeighbor,
      noZlowerNeighbor, noZupperNeighbor);
    // if the particle is exiting, put it in the appropriate send bin;
    // this could be done at conclusion of push after particles are
    // converted to AoS format in order to overlap communication
    // with computation.
    bool was_sent = send_pcl_to_appropriate_buffer(pcl,
      isBoundaryProcess,
      hasXlowerNeighbor, hasXupperNeighbor,
      hasYlowerNeighbor, hasYupperNeighbor,
      hasZlowerNeighbor, hasZupperNeighbor);
    if(was_sent)
    {
      _pcls.delete_element(np_current);
    }
    else
    {
      np_current++;
    }
  }
  assert_eq(_pcls.size(),np_current);

  // receive and redistribute particles once for
  // each dimension of space without doing an
  // all-reduce to check if any particles are
  // actually being communicated.
  int num_pcls_resent;
  for(int i=0;i<3;i++)
  {
    flush_send();
    num_pcls_resent = handle_incoming_particles();
  }

  // continue receiving and resending incoming particles until
  // global all-reduce of num_pcls_resent is zero, indicating
  // that there are no more particles to be received.
  //
  long long total_num_pcls_resent = mpi_global_sum(num_pcls_resent);
  while(total_num_pcls_resent)
  {
    flush_send();
    num_pcls_resent = handle_incoming_particles();
    total_num_pcls_resent = mpi_global_sum(num_pcls_resent);
  }
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
  VirtualTopology3D *ptVCT = vct;
  cout << endl;
  cout << "Number of Particles: " << _pcls.size() << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << get_species_num() << endl;
  for (int i = 0; i < _pcls.size(); i++)
  {
    SpeciesParticle& pcl = _pcls[i];
    const double x = pcl.get_x();
    const double y = pcl.get_y();
    const double z = pcl.get_z();
    cout << "Particle #" << i << ": x=" << x << " y=" << y << " z=" << z << " u=" << u << " v=" << v << " w=" << w << endl;
  }
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp()  const
{
  VirtualTopology3D *ptVCT = vct;
  cout << endl;
  cout << "Number of Particles of species " << get_species_num() << ": " << getNOP() << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
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
      updateParticlesToSoA();
      break;
    default:
      unsupported_value_error(particleType);
  }
}

// need to sort and communicate particles after each iteration
void Particles3Dcomm::sort_particles_serial_AoS()
{
  convertParticlesToAoS();

  Larray<SpeciesParticle>& pcls = fetch_pcls();
  Larray<SpeciesParticle>& pclstmp = fetch_pclstmp();
  pclstmp.reserve(pcls.size());
  {
    numpcls_in_bucket->setall(0);
    // iterate through particles and count where they will go
    for (int pidx = 0; pidx < pcls.size(); pidx++)
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
    assert_eq(accpcls,nop);

    numpcls_in_bucket_now->setall(0);
    // put the particles where they are supposed to go
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
      pclstmp[outpidx] = pcl;
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
  u.resize(nop);
  v.resize(nop);
  w.resize(nop);
  q.resize(nop);
  x.resize(nop);
  y.resize(nop);
  z.resize(nop);
  t.resize(nop);
 #ifndef __MIC__
  #pragma omp for
  for(int pidx=0; nop; pidx++)
  {
    const SpeciesParticle& pcl = pcls[pidx];
    u[pidx] = pcl.get_u(0);
    v[pidx] = pcl.get_u(1);
    w[pidx] = pcl.get_u(2);
    q[pidx] = pcl.get_q();
    x[pidx] = pcl.get_x(0);
    y[pidx] = pcl.get_x(1);
    z[pidx] = pcl.get_x(2);
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
  particleType = ParticleType::both;
}

void Particles3Dcomm::copyParticlesToAoS()
{
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_AOS);
  _pcls.resize(u.size());
  dprintf("copying to array of structs");
 #ifndef __MIC__
  // use a simple stride-8 gather
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx++)
  {
    pcls[pidx].set(
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
  particleType = ParticleType::both;
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

// defines SoA to be the authority
// (conceptually releasing any write-lock)
//
bool Particles3Dcomm::particlesAreSoA()
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

