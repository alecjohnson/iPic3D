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
#include "VCtopology3D.h"
#include "Collective.h"
#include "Alloc.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid3DCU.h"
#include "MPIdata.h"
#include "ompdefs.h"
#include "ipic_math.h"
#include "ipic_defs.h"
#include "mic_basics.h"
#include "parallel.h"

#include "Particle.h"
#include "Particles3Dcomm.h"
#include "Parameters.h"

#include "ipic_hdf5.h"
#include "Restart3D.h"
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

static bool print_pcl_comm_counts = false;

static void print_pcl(SpeciesParticle& pcl, int is)
{
  dprintf("--- pcl spec %d ---", is);
  dprintf("u = %+6.4f", pcl.get_u());
  dprintf("v = %+6.4f", pcl.get_v());
  dprintf("w = %+6.4f", pcl.get_w());
  dprintf("q = %+6.4f", pcl.get_q());
  dprintf("x = %+6.4f", pcl.get_x());
  dprintf("y = %+6.4f", pcl.get_y());
  dprintf("z = %+6.4f", pcl.get_z());
  dprintf("t = %5.0f", pcl.get_t());
}

static void print_pcls(vector_SpeciesParticle& pcls, int start, int is)
{
  for(int pidx=start; pidx<pcls.size();pidx++)
  {
    dprintf("--- particle %d.%d ---", is,pidx);
    dprintf("u[%d] = %+6.4f", pidx, pcls[pidx].get_u());
    dprintf("v[%d] = %+6.4f", pidx, pcls[pidx].get_v());
    dprintf("w[%d] = %+6.4f", pidx, pcls[pidx].get_w());
    dprintf("q[%d] = %+6.4f", pidx, pcls[pidx].get_q());
    dprintf("x[%d] = %+6.4f", pidx, pcls[pidx].get_x());
    dprintf("y[%d] = %+6.4f", pidx, pcls[pidx].get_y());
    dprintf("z[%d] = %+6.4f", pidx, pcls[pidx].get_z());
    dprintf("t[%d] = %5.0f", pidx, pcls[pidx].get_t());
  }
}
void print_pcls(vector_SpeciesParticle& pcls, int is, longid* id_list, int num_ids)
{
  dprintf("=== species %d, with %d pcls ===", is, pcls.size());
  for(int pidx=0; pidx<pcls.size();pidx++)
  for(int i=0;i<num_ids;i++)
  if(pcls[pidx].get_ID()==id_list[i])
  {
    dprintf("--- particle %d.%d ---", is,pidx);
    dprintf("u[%d] = %+6.4f", pidx, pcls[pidx].get_u());
    dprintf("v[%d] = %+6.4f", pidx, pcls[pidx].get_v());
    dprintf("w[%d] = %+6.4f", pidx, pcls[pidx].get_w());
    dprintf("q[%d] = %+6.4f", pidx, pcls[pidx].get_q());
    dprintf("x[%d] = %+6.4f", pidx, pcls[pidx].get_x());
    dprintf("y[%d] = %+6.4f", pidx, pcls[pidx].get_y());
    dprintf("z[%d] = %+6.4f", pidx, pcls[pidx].get_z());
    dprintf("t[%d] = %5.0f", pidx, pcls[pidx].get_t());
  }
}

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
  const Setting& setting_)
 :
  is(species_number),
  setting(setting_),
  vct(&setting.vct()),
  grid(&setting.grid()),
  pclIDgenerator(),
  particleType(ParticleType::AoS)
{
  // communicators for particles
  //
  MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
  //
  // define connections
  //using namespace Direction;
  //
  // reserve one communicator for neighbor in every direction
  sendPclList.reserve(NUM_COMM_NEIGHBORS);
  recvPclList.reserve(NUM_COMM_NEIGHBORS);
  //const int *mycoords = vct->getCoordinates();
  // for each direction index:
  // * recvPclList Connection is populated from sources and directions:
  //   - direction tag is uncorrected displacement
  //     from potential source to destination
  //   - source is results of pinning potential source.
  // * sendPclList Connection is populated from sources and directions:
  //   - direction tag is pin-reflection of uncorrected displacement
  //     from source to potential destination (thus yielding displacement
  //     from potential source to destination
  //   - destination is result of pinning potential destination.
  // * adding direction to mycoords gives potential source
  // 
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
  {
    // for each direction determine the destination process...
    int dest_rank = vct->get_valid_neighbor_rank(di);
    // ...and an appropriate direction tag.
    int send_tag = vct->get_pinned_direction_tag(di);
    // the direction tag encodes the direction that the data moves
    sendPclList[di].init(Connection(dest_rank, send_tag, mpi_comm));
  }
  // di is the direction that the data is coming from
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
  {
    // from_di is the direction that the data is moving
    int from_di = vct->get_negated_direction_tag(di);
    // for each direction determine the destination process...
    int source_rank = vct->get_valid_neighbor_rank(from_di);
    // the direction tag encodes the direction that the data moves
    recvPclList[di].init(Connection(source_rank, di, mpi_comm));
  }

  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
  {
    recvPclList[di].post_recvs();
  }

  // info from collectiveIO
  //
  const Collective* col = &setting.col();
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
  beg_pos[0] = grid->getXstart();
  beg_pos[1] = grid->getYstart();
  beg_pos[2] = grid->getZstart();
  end_pos[0] = grid->getXend();
  end_pos[1] = grid->getYend();
  end_pos[2] = grid->getZend();
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
  if (restart != 0)
  {
  #ifdef NO_HDF5
    eprintf("restart is supported only if compiling with HDF5");
  #else
    int species_number = get_species_num();
    // prepare arrays to receive particles
    particleType = ParticleType::SoA;
    read_particles_restart(col, vct, species_number,
      u, v, w, q, x, y, z, t);
    convertParticlesToAoS();
  #endif
  }

  // velocity capping
  //
  const double domain_crossing_vel[3] = {
    col->getLx()/col->getDt(),
    col->getLy()/col->getDt(),
    col->getLz()/col->getDt(),
  };
  double subdomain_crossing_vel[3] = {
    domain_crossing_vel[0]/col->getXLEN(),
    domain_crossing_vel[1]/col->getYLEN(),
    domain_crossing_vel[2]/col->getZLEN(),
  };
  double cell_crossing_vel[3] = {
    col->getDx()/col->getDt(),
    col->getDy()/col->getDt(),
    col->getDz()/col->getDt(),
  };
  switch(Parameters::vel_cap_scale_ref())
  {
    using namespace Parameters;
    default:
      unsupported_value_error(vel_cap_scale_ref());
    case null: // default initialization of maxvel
    case fraction_of_subdomain:
    {
      maxvel[0] = subdomain_crossing_vel[0];
      maxvel[1] = subdomain_crossing_vel[1];
      maxvel[2] = subdomain_crossing_vel[2];
      break;
    }
    case fraction_of_domain:
    {
      maxvel[0] = domain_crossing_vel[0];
      maxvel[1] = domain_crossing_vel[1];
      maxvel[2] = domain_crossing_vel[2];
      break;
    }
    case absolute_velocity:
      maxvel[0] = maxvel[1] = maxvel[2] = 1.;
      break;
    case fraction_of_light_speed:
      maxvel[0] = maxvel[1] = maxvel[2] = col->getC();
      break;
  }
  maxvel[0] *= Parameters::vel_cap_scaled_value();
  maxvel[1] *= Parameters::vel_cap_scaled_value();
  maxvel[2] *= Parameters::vel_cap_scaled_value();
  minvel[0] = -maxvel[0];
  minvel[1] = -maxvel[1];
  minvel[2] = -maxvel[2];
  // fraction of domain that particle can move
  double subdomain_fractional_vel[3];
  // maximum number of communications needed to communicate all particles.
  double max_num_pcl_comms_double = 0.9; // 1
  for(int i=0;i<3;i++)
  {
    subdomain_fractional_vel[i] = maxvel[i]/subdomain_crossing_vel[i];
    max_num_pcl_comms_double
      = std::max(max_num_pcl_comms_double,subdomain_fractional_vel[i]);
  }
  max_num_pcl_comms = ceil(max_num_pcl_comms_double);
  if(!Parameters::vel_cap_scale_ref())
    assert_eq(max_num_pcl_comms,1);
  // show velocity cap that will be applied
  if(is_output_thread())
  {
    printf("--- particle mover parameters ---\n");
    printf("  c=%g\n", col->getC());
    printf("  species %d velocity cap: umax=%g,vmax=%g,wmax=%g\n",
      is, maxvel[0],maxvel[1],maxvel[2]);
    printf("  species %d thermal vel: uth=%g,vth=%g,wth=%g\n",
      is, col->getUth(is), col->getVth(is), col->getWth(is));
    printf("  cell crossing velocities: [%g,%g,%g]\n",
      cell_crossing_vel[0],
      cell_crossing_vel[1],
      cell_crossing_vel[2]);
    printf("  subdomain crossing velocities: [%g,%g,%g]\n",
      subdomain_crossing_vel[0],
      subdomain_crossing_vel[1],
      subdomain_crossing_vel[2]);
    printf("  particles should be communicated in at most %d=ceil(%g)"
      " communications\n", max_num_pcl_comms, max_num_pcl_comms_double);
  }
}

bool Particles3Dcomm::print_pcl_comm_counts()const
{
  if(is==0 && ::print_pcl_comm_counts)
    return true;
  return false;
}
// pad capacities so that aligned vectorization
// does not result in an array overrun.
//
// this should usually be cheap (a no-op)
//
void Particles3Dcomm::pad_capacities()
{
 #pragma omp single
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
}

void Particles3Dcomm::resize_AoS(int nop)
{
 #pragma omp master
 {
  const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
  _pcls.reserve(padded_nop);
  _pcls.resize(nop);
 }
}

void Particles3Dcomm::resize_SoA(int nop)
{
 #pragma omp master
 {
  //
  // allocate space for particles including padding
  //
  const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
  if(is_output_thread()) dprintf("allocating to hold %d", padded_nop);
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
  if(is_output_thread()) dprintf("done resizing to hold %d", nop);
 }
}

// returns true if particle was sent
//
// should vectorize this by comparing position vectors
//
inline bool Particles3Dcomm::send_pcl_to_appropriate_buffer(
  SpeciesParticle& pcl, int*send_count)
{
  const double *u = &pcl.fetch_u();
  const double *beg_pos = grid->get_beg_pos();
  const double *end_pos = grid->get_end_pos();

  ASSUME_ALIGNED(u);
  ASSUME_ALIGNED(beg_pos);
  ASSUME_ALIGNED(end_pos);
  // which direction is the particle going?
  double x[3] ALLOC_ALIGNED;
  int dirs[3] ALLOC_ALIGNED;
  for(int i=0;i<3;i++)
  {
    x[i] = u[4+i];
    // this way would avoid branching but does not cap direction
    //int dirs[i] = floor((x[i]-beg_pos[i])*widthinv[i])
    if(x[i]<beg_pos[i]) dirs[i]=-1;
    else if(x[i]>end_pos[i]) dirs[i]=1;
    else dirs[i]=0;
  }
  // if going nowhere then return
  bool send_pcl = dirs[0]||dirs[1]||dirs[2];
  if(!send_pcl)
    return false;

  // put particle in appropriate communication buffer if exiting
  int dir_idx = vct->get_index_for_direction(dirs);
  // send particle
  sendPclList[dir_idx].send(pcl);
  send_count[dir_idx]++;
  return true;;
}

// flush sending particles.
//
void Particles3Dcomm::flush_send()
{
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
    sendPclList[di].send_complete();
}

void Particles3Dcomm::apply_periodic_BC_global(
  vector_SpeciesParticle& pcl_list, int pstart)
{
  double Ls[3] ALLOC_ALIGNED;
  Ls[0] = setting.col().getLx();
  Ls[1] = setting.col().getLy();
  Ls[2] = setting.col().getLz();
  double Linvs[3] ALLOC_ALIGNED;
  for(int i=0;i<3;i++) Linvs[i] = 1/Ls[i];

  // apply shift to all periodic directions
  for(int pidx=pstart;pidx<pcl_list.size();pidx++)
  for(int i=0;i<3;i++)
  {
    SpeciesParticle& pcl = pcl_list[pidx];
    double *x = &pcl.fetch_x();
    if(vct->getPeriods(i))
    {
      x[i] = modulo(x[i], Ls[i], Linvs[i]);
    }
  }
}

// routines for sorting list of particles
//
// sort_pcls: macro to put all particles that satisfy a condition
// at the end of an array of given size.  It has been written so
// that it could be replaced with a generic routine that takes
// the "condition" method as a callback function (and if the
// optimizer is good, the performance in this case would actually
// be just as good, so maybe we should just make this change).
//
// pcls: SpeciesParticle* or vector_SpeciesParticle
// size_in: number of elements in pcls list
// start_in: starting index of list to be sorted
// start_out: returns starting index of particles
//    for which the condition is true.
// condition: a (probably inline) function of SpeciesParticle
//   that returns true if the particle should go at the end of the list
//
#define sort_pcls(pcls, start_in, start_out, condition) \
{ \
  int start = (start_in); \
  assert(0<=start); \
  start_out = pcls.size(); \
  /* pidx traverses the array */ \
  for(int pidx=pcls.size()-1;pidx>=start;pidx--) \
  { \
    assert(pidx<start_out); \
    /* if condition is true, put the particle at the end of the list */ \
    if(condition(pcls[pidx])) \
    { \
      --start_out; \
      SpeciesParticle tmp_pcl = pcls[pidx]; \
      pcls[pidx] = pcls[start_out]; \
      pcls[start_out] = tmp_pcl; \
    } \
  } \
}
//  /* pidx increases and start_out decreases until they meet */ \
//  for(int pidx=start;pidx<start_out;) \
//  { \
//    /* if condition is true, put the particle at the end of the list */ \
//    if(condition(pcls[pidx])) \
//    { \
//      --start_out; \
//      SpeciesParticle tmp_pcl = pcls[pidx]; \
//      pcls[pidx] = pcls[start_out]; \
//      pcls[start_out] = tmp_pcl; \
//    } \
//    else \
//    { \
//      pidx++; \
//    } \
//  } \

// condition methods to use in sorting particles
inline bool Particles3Dcomm::test_outside_domain(const SpeciesParticle& pcl)const
{
  // This could be vectorized
  bool is_outside_domain=(
       pcl.get_x() < 0. || pcl.get_y() < 0. || pcl.get_z() < 0.
    || pcl.get_x() > Lx || pcl.get_y() > Ly || pcl.get_z() > Lz );
  return is_outside_domain;
}
inline bool Particles3Dcomm::test_outside_nonperiodic_domain(const SpeciesParticle& pcl)const
{
  // This could be vectorized
  bool is_outside_nonperiodic_domain =
     (!vct->getPeriods(0) && (pcl.get_x() < 0. || pcl.get_x() > Lx)) ||
     (!vct->getPeriods(1) && (pcl.get_y() < 0. || pcl.get_y() > Ly)) ||
     (!vct->getPeriods(2) && (pcl.get_z() < 0. || pcl.get_z() > Lz));
  return is_outside_nonperiodic_domain;
}

// apply user-supplied boundary conditions
//
// A better implementation would sort particles
// into 27 categories based on their direction index
// and then call a version of apply_BCs
// to each segment.
void Particles3Dcomm::apply_nonperiodic_BCs_global(
  vector_SpeciesParticle& pcl_list, int pstart)
{
  int lstart;
  int lsize;
  if(!vct->getPeriods(0))
  {
    // separate out particles that need Xleft boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Xleft_of_domain);
    // apply boundary conditions
    apply_Xleft_BC(pcl_list, lstart);
    // separate out particles that need Xrght boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Xrght_of_domain);
    // apply boundary conditions
    apply_Xrght_BC(pcl_list, lstart);
  }
  if(!vct->getPeriods(1))
  {
    // separate out particles that need Yleft boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Yleft_of_domain);
    // apply boundary conditions
    apply_Yleft_BC(pcl_list, lstart);
    // separate out particles that need Yrght boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Yrght_of_domain);
    // apply boundary conditions
    apply_Yrght_BC(pcl_list, lstart);
  }
  if(!vct->getPeriods(2))
  {
    // separate out particles that need Zleft boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Zleft_of_domain);
    // apply boundary conditions
    apply_Zleft_BC(pcl_list, lstart);
    // separate out particles that need Zrght boundary conditions applied
    sort_pcls(pcl_list, pstart, lstart, test_Zrght_of_domain);
    // apply boundary conditions
    apply_Zrght_BC(pcl_list, lstart);
  }
}

bool Particles3Dcomm::test_pcls_are_in_nonperiodic_domain(const vector_SpeciesParticle& pcls)const
{
  const int size = pcls.size();
  for(int pidx=0;pidx<size;pidx++)
  {
    const SpeciesParticle& pcl = pcls[pidx];
    // should vectorize these comparisons
    bool not_in_domain = test_outside_nonperiodic_domain(pcl);
    if(__builtin_expect(not_in_domain, false)) return false;
  }
  return true; // all pcls are in domain
}
bool Particles3Dcomm::test_pcls_are_in_domain(const vector_SpeciesParticle& pcls)const
{
  const int size = pcls.size();
  for(int pidx=0;pidx<size;pidx++)
  {
    const SpeciesParticle& pcl = pcls[pidx];
    // should vectorize these comparisons
    bool not_in_domain = test_outside_domain(pcl);
    if(__builtin_expect(not_in_domain, false)) return false;
  }
  return true; // all pcls are in domain
}
bool Particles3Dcomm::test_all_pcls_are_in_subdomain()
{
  const int size = _pcls.size();
  for(int pidx=0;pidx<size;pidx++)
  {
    SpeciesParticle& pcl = _pcls[pidx];
    // should vectorize these comparisons
    bool not_in_subdomain = 
         pcl.get_x() < beg_pos[0] || pcl.get_y() < beg_pos[1] || pcl.get_z() < beg_pos[2]
      || pcl.get_x() > end_pos[0]   || pcl.get_y() > end_pos[1]   || pcl.get_z() > end_pos[1];
    if(__builtin_expect(not_in_subdomain, false)) return false;
  }
  return true; // all pcls are in processor subdomain
}

// If do_apply_periodic_BC_global is false, then we may need to
// communicate particles as many as 2*(XLEN+YLEN+ZLEN) times
// after the call to handle_received_particles(true), because it
// is conceivable that a particle will be communicated the full
// length of a periodic domain, then have its position remapped
// to the proper position, and then traverse the full domain
// again.
//
// If do_apply_periodic_BC_global is true, then we can
// guarantee that all particles will be communicated within
// at most (XLEN+YLEN+ZLEN) communications, but on the other
// hand, a particle that would have been communicated in only
// two iterations by being wrapped around a periodic boundary
// will instead be communicated almost the full length of that
// dimension.
static bool do_apply_periodic_BC_global = false;
// apply boundary conditions to a list of particles globally
// (i.e. without regard to their current location in memory)
void Particles3Dcomm::apply_BCs_globally(vector_SpeciesParticle& pcl_list)
{
  // apply boundary conditions to every
  // particle until every particle lies in the domain.

  // put the particles outside the domain at the end of the list
  //
  // index of first particle that is unfinished
  int pstart = 0;
  // sort particles outside of the domain to the end of the list
  sort_pcls(pcl_list, 0, pstart, test_outside_domain);

  for(int i=0; pstart < pcl_list.size(); i++)
  {
    if(do_apply_periodic_BC_global)
    {
      apply_periodic_BC_global(pcl_list, pstart);
      // apply user-supplied boundary conditions
      apply_nonperiodic_BCs_global(pcl_list, pstart);
      // put particles outside of the domain at the end of the list
      sort_pcls(pcl_list, pstart, pstart, test_outside_domain);
    }
    else
    {
      apply_nonperiodic_BCs_global(pcl_list, pstart);
      sort_pcls(pcl_list, pstart, pstart, test_outside_nonperiodic_domain);
    }

    // if this fails, something has surely gone wrong
    // (e.g. we have a runaway particle).
    if(i>=100)
    {
      print_pcls(pcl_list,pstart, get_species_num());
      dprint(pstart);
      dprint(pcl_list.size());
      eprintf("something went wrong.")
    }
  }
  if(do_apply_periodic_BC_global)
  {
    assert(test_pcls_are_in_domain(pcl_list));
  }
  else
  {
    assert(test_pcls_are_in_nonperiodic_domain(pcl_list));
  }
}

// dir_idx is the direction that the data is moving,
// so the direction that the data is coming from is -dir.
//
// Note that the documentation in the input file says that
// boundary conditions are simply ignored in the periodic
// case, so no check is made.
//
void Particles3Dcomm::apply_BCs_locally(vector_SpeciesParticle& pcl_list,
  int dir_idx)
{
  // determine direction of particle transfer
  int dirs[3];
  vct->set_direction_for_index(dirs, dir_idx);
  int src_dirs[3]={ -dirs[0], -dirs[1], -dirs[2] };
  // determine if boundary conditions need to be applied
  // and if so for which directions
  int applyBCdirs[3]={0,0,0};
  for(int i=0;i<3;i++)
  {
    if((vct->noLeftNeighbor(i) && src_dirs[i]<0) ||
       (vct->noRghtNeighbor(i) && src_dirs[i]>0))
    {
      applyBCdirs[i] = src_dirs[i];
    }
  }
  const bool do_apply_BCs = applyBCdirs[0] || applyBCdirs[1] || applyBCdirs[2];

  // determine if periodicity shift is needed
  // and if so for which directions
  bool wrap_dir[3]={false,false,false};
  for(int i=0;i<3;i++)
  {
    if((vct->isPeriodicLower(i) && src_dirs[i]<0) ||
       (vct->isPeriodicUpper(i) && src_dirs[i]>0))
    {
      wrap_dir[i]=true;
    }
  }
  const bool apply_shift = wrap_dir[0] || wrap_dir[1] || wrap_dir[2];

  // if appropriate apply periodicity shift for this block
  //
  // could change from modulo to simple shift.
  //
  double Ls[3] ALLOC_ALIGNED;
  Ls[0] = setting.col().getLx();
  Ls[1] = setting.col().getLy();
  Ls[2] = setting.col().getLz();
  double Linvs[3] ALLOC_ALIGNED;
  for(int i=0;i<3;i++) Linvs[i] = 1/Ls[i];

  if(apply_shift)
  {
    for(int pidx=0;pidx<pcl_list.size();pidx++)
    {
      SpeciesParticle& pcl = pcl_list[pidx];
      double *xptr = &pcl.fetch_x();
      for(int i=0;i<3;i++)
      {
        if(wrap_dir[i])
        {
          xptr[i] = modulo(xptr[i], Ls[i], Linvs[i]);
          // The following would be more efficient and would
          // be good enough if we cap particle motion, but
          // when solving 2D problems with only one layer of
          // mesh cells in an ignorable Z direction the user
          // must be sure to choose Lz sufficiently large to
          // ensure that the fastest particles do not cross
          // the Z dimension of the domain.
          //xptr[i] += dirs[i]*Ls[i];
        }
      }
    }
  }
  // if appropriate then apply boundary conditions to this block
  if(do_apply_BCs)
  {
    apply_BCs(pcl_list, applyBCdirs);
    // this was an open boundary conditions bug:
    //int size = pcl_list.size();
    //pcl_list.resize(size);
  }
}

namespace PclCommMode
{
  enum Enum
  {
    do_apply_BCs_globally=1,
    print_sent_pcls=2,
  };
}
// receive, sort, and, as appropriate, resend incoming particles
//
// assumes that flush_send() has been called
//
// returns number of particles that were resent
//
int Particles3Dcomm::handle_received_particles(int pclCommMode)
{
  using namespace PclCommMode;

  // reinitialize communication
  //
  int recv_count[NUM_COMM_NEIGHBORS];
  int send_count[NUM_COMM_NEIGHBORS];
  MPI_Request recv_requests[NUM_COMM_NEIGHBORS];
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
  {
    // we expect to receive at least one block from every
    // communicator, so make sure that all receive buffers are
    // clear and waiting
    recvPclList[di].recv_start();
    // make sure that current block in each sender is ready for sending
    sendPclList[di].send_start();
    recv_count[di]=0;
    send_count[di]=0;
    recv_requests[di]=recvPclList[di].get_curr_request();
    // confirm that no communicators are finished
    assert(!recvPclList[di].comm_finished());
  }

  int num_pcls_recved = 0;
  int num_pcls_resent = 0;
  // we expect communication from all communicators.
  int num_unfinished_communicators = NUM_COMM_NEIGHBORS;
  // while there are still incoming particles
  // put them in the appropriate buffer.
  while(num_unfinished_communicators)
  {
    // wait until one of the receiving buffers finishes receiving a message
    //
    int recv_index;
    MPI_Status recv_status;
    pcls_MPI_Waitany(NUM_COMM_NEIGHBORS,recv_requests,&recv_index,&recv_status);
    //
    // check that the received message data makes sense
    //
    if(recv_index==MPI_UNDEFINED)
      eprintf("recv_requests contains no active handles");
    assert_ge(recv_index,0);
    assert_lt(recv_index,NUM_COMM_NEIGHBORS);
    //
    // grab the received block of particles and process it
    //
    BlockCommunicator<SpeciesParticle>& recvComm = recvPclList[recv_index];
    // this act of fetching determines whether communication is finished
    // and number of particles received.
    Block<SpeciesParticle>& recv_block
      = recvComm.fetch_received_block(recv_status);
    // track finished communicators
    if(recvComm.comm_finished())
      num_unfinished_communicators--;
    // track number of particles received
    recv_count[recv_index]+=recv_block.size();
    num_pcls_recved += recv_block.size();
    // grab the list of particles to process
    vector_SpeciesParticle& pcl_list = recv_block.fetch_block();

    // apply boundary conditions to the incoming particles
    if(pclCommMode&do_apply_BCs_globally)
    {
      apply_BCs_globally(pcl_list);
    }
    else
    {
      // recv_index is the index of the direction we are receiving from
      apply_BCs_locally(pcl_list, recv_index);
    }

    // process each particle in the received block.
    {
      for(int pidx=0;pidx<recv_block.size();pidx++)
      {
        SpeciesParticle& pcl = recv_block[pidx];
        bool was_sent = send_pcl_to_appropriate_buffer(pcl, send_count);

        if(__builtin_expect(was_sent,false))
        {
          num_pcls_resent++;
          if(pclCommMode&print_sent_pcls)
          {
            print_pcl(pcl,get_species_num());
          }
        }
        else
        {
          // The particle belongs here, so put it in the
          // appropriate place. For now, all particles are in a
          // single list, so we append to the list.
          _pcls.push_back(pcl);
        }
      }
    }
    // release the block and update the receive request
    recvComm.release_received_block();
    recv_requests[recv_index] = recvComm.get_curr_request();
  }
  // confirm that all communicators are finished
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
    assert(recvPclList[di].comm_finished());

  // report the number of particles received e.g. per direction
  if(print_pcl_comm_counts())
  {
    dprintf("spec %d recved_count: %d", is, num_pcls_recved);
    dprintf("spec %d resent_count: %d", is, num_pcls_resent);
    for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
    {
      int dirs[3];
      vct->set_direction_for_index(dirs, di);
      dprintf("spec %d recv in dir [%d,%d,%d]: %d", is,
        dirs[0], dirs[1], dirs[2], recv_count[di]);
    }
  }
  // return the number of particles that were resent

  return num_pcls_resent;
}

static long long mpi_global_sum(int in)
{
  long long total;
  long long long_in = in;
  //dprintf("calling MPI_Allreduce(%d,&total,1, ...)", long_in);
  timeTasks_set_task(TimeTasks::PCLS_MPI_ALLREDUCE);
  MPI_Allreduce(&long_in, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  return total;
}

// this methods should be virtual
// so that the user can override boundary conditions.
//
// apply BCs to all particles in pcls starting with
// index start.
//
// It is expected that dirs=[xoffset, yoffset, zoffset],
// where xoffset is -1 if particles are to the left of the
// process subdomain, 1 if to the right of the process
// subdomain, and 0 if in the process subdomain.
// (The same situation is assumed to hold for all particles
// in the list).
//
void Particles3Dcomm::apply_BCs(vector_SpeciesParticle& pcls,
  const int dirs[3], int start)
{
  // if the particles are in the domain then there is nothing to do
  if(!(dirs[0] || dirs[1] || dirs[2]))
  {
    // should not get here
    assert(false);
    return;
  }

  // if the particles exited an open boundary then delete them
  const bool exiting_domain = 
    (bcPfaceXleft == BCparticles::EXIT && dirs[0] < 0) ||
    (bcPfaceXright== BCparticles::EXIT && dirs[0] > 0) ||
    (bcPfaceYleft == BCparticles::EXIT && dirs[1] < 0) ||
    (bcPfaceYright== BCparticles::EXIT && dirs[1] > 0) ||
    (bcPfaceZleft == BCparticles::EXIT && dirs[2] < 0) ||
    (bcPfaceZright== BCparticles::EXIT && dirs[2] > 0);

  if(exiting_domain)
  {
    assert(false); // this code has not yet been tested.
    pcls.resize(start);
    return;
  }

  // parameter digestion
  //
  int size = pcls.size();
  assert_le(0,start);
  // reemission reverses direction of velocity
  const int newdirs[3] = {-dirs[0],-dirs[1],-dirs[2]};

  // if the particles crossed a reemission boundary then reemit them
  const bool do_reemit = 
    (bcPfaceXleft == BCparticles::REEMISSION && dirs[0] < 0) ||
    (bcPfaceXright== BCparticles::REEMISSION && dirs[0] > 0) ||
    (bcPfaceYleft == BCparticles::REEMISSION && dirs[1] < 0) ||
    (bcPfaceYright== BCparticles::REEMISSION && dirs[1] > 0) ||
    (bcPfaceZleft == BCparticles::REEMISSION && dirs[2] < 0) ||
    (bcPfaceZright== BCparticles::REEMISSION && dirs[2] > 0);

  // Algorithm:
  //
  // 1. apply reemission if a reemission boundary is involved
  //
  if(do_reemit)
  {
    assert(false); // this code has not yet been tested.
    for(int p=start;p<size;p++)
    {
      SpeciesParticle& pcl = pcls[p];
      double *uptr = &pcl.fetch_u();
      sample_maxwellian(uptr[0],uptr[1],uptr[2], uth,vth,wth);
    }
  }
  //
  // 2. reflect positions and velocities appropriately
  //
  const Collective& col = setting.col();
  const double Ls[3]={col.getLx(),col.getLy(),col.getLz()};
  double signs[3] ALLOC_ALIGNED;
  double wall2[3] ALLOC_ALIGNED;
  for(int i=0;i<3;i++)
  {
    // reflect unless in middle
    signs[i] = dirs[i] ? -1 : 1;
    wall2[i] = dirs[i]>0 ? 2*Ls[i] : 0;
  }
  // the order of these two loops should be swapped
  // unless the compiler decides to vectorize.
  for(int p=start;p<size;p++)
  for(int i=0;i<3;i++)
  {
    SpeciesParticle& pcl = pcls[p];
    double *xptr = &pcl.fetch_x();
    double *uptr = &pcl.fetch_u();
    if(dirs[i])
    {
      // this formula reflects off wall if dir is nonzero
      // and would do nothing is dir is zero
      xptr[i] = wall2[i] + signs[i]*xptr[i];
      uptr[i] = newdirs[i]*fabs(uptr[i]);
      // old way (but then reemission must
      // handle reflection separately).
      // uptr[i] *= -1;
    }
  }
}

// these methods should be made virtual
// so that the user can override boundary conditions.
//
void Particles3Dcomm::apply_Xleft_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={-1,0,0};
  apply_BCs(pcls, dirs, start);
}
void Particles3Dcomm::apply_Yleft_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={0,-1,0};
  apply_BCs(pcls, dirs, start);
}
void Particles3Dcomm::apply_Zleft_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={0,0,-1};
  apply_BCs(pcls, dirs, start);
}
void Particles3Dcomm::apply_Xrght_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={1,0,0};
  apply_BCs(pcls, dirs, start);
}
void Particles3Dcomm::apply_Yrght_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={0,1,0};
  apply_BCs(pcls, dirs, start);
}
void Particles3Dcomm::apply_Zrght_BC(vector_SpeciesParticle& pcls, int start)
{
  const int dirs[3]={0,0,1};
  apply_BCs(pcls, dirs, start);
}

// return number of particles sent
int Particles3Dcomm::separate_and_send_particles()
{
  // find particle 0 and print it
  if(false) {
    const int num_ids = 1;
    longid id_list[num_ids] = {0};
    print_pcls(_pcls,get_species_num(),id_list, num_ids);
  }
  timeTasks_set_communicating(); // communicating until end of scope

  convertParticlesToAoS();

  //
  int send_count[NUM_COMM_NEIGHBORS];
  for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
  {
    send_count[di]=0;
    // activate receiving
    recvPclList[di].recv_start();
    // make sure that current block is ready for sending
    sendPclList[di].send_start();
  }

  //int send_count[6]={0,0,0,0,0,0};
  const int num_pcls_initially = _pcls.size();
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

  // flush all send buffers
  flush_send();

  // check and report statistics about number of particles sent
  //
  assert_eq(_pcls.size(),np_current);
  const int num_pcls_unsent = getNOP();
  const int num_pcls_sent = num_pcls_initially - num_pcls_unsent;
  if(print_pcl_comm_counts())
  {
    for(int di=0;di<NUM_COMM_NEIGHBORS;di++)
    {
      int dirs[3];
      vct->set_direction_for_index(dirs, di);
      dprintf("spec %d send_count in dir [%d,%d,%d]: %d",get_species_num(),
        dirs[0], dirs[1], dirs[2], send_count[di]);
    }
    dprintf("spec %d total_sent = %d",is,num_pcls_sent);
  }
  return num_pcls_sent;
}

// communicate particles and apply boundary conditions
// until every particle is in the process of its subdomain.
//
// At the end of this method, the position of every particle in
// this process must lie in this process's proper subdomain,
// for two reasons: (1) sumMoments assumes that all particles
// lie in the proper subdomain of the process, and if this
// assumption is violated memory corruption can result, and
// (2) the extrapolation algorithm used by the mover is unstable,
// which could cause the velocity of particles not in the proper
// subdomain to blow up.
//
// The algorithm proceeds in three steps:
// (1) for min_num_iterations, receive particles from neighbor
// processes, apply boundary conditions locally (i.e. only if
// this is a boundary process), and resend any particles that
// still do not belong in this subdomain,
// (2) apply boundary conditions globally (i.e. independent of
// whether this is a boundary process), and resend any particles
// not in this subdomain, and
// (3) communicate particles (as many as XLEN+YLEN+ZLEN times)
// until every particle is in its appropriate subdomain.
//
// This communication algorithm is perhaps more sophisticated
// than is justified by the mover, which does not properly
// resolve fast-moving particles.
//
// min_num_iterations: number of iterations that this process
//   applies boundary conditions locally.
//   This method then applies boundary conditions globally to each
//   incoming particle until the particle resides in the domain.
//   forces
//   
void Particles3Dcomm::receive_particles_until_done()
{
  timeTasks_set_communicating(); // communicating until end of scope
  assert_gt(max_num_pcl_comms,0);
  //
  // enforcing cap on particle velocities determines
  // a cap on the number of iterations needed.
  // only one iteration is needed if particles
  // can move at most one subdomain per time step.
  //
  int num_pcls_sent;
  //
  // we prefer to do a fixed number of local
  // communications rather than performing an
  // all-reduce on the number of particles sent.
  //
  int num_recommunications = max_num_pcl_comms-1;
  for(int i=0;i<num_recommunications;i++)
  {
    num_pcls_sent = handle_received_particles();
    flush_send(); // flush sending of particles
  }
  // handle the received particles the final (and perhaps only) time
  num_pcls_sent = handle_received_particles();

  if(Parameters::vel_cap_method()==Parameters::box_capped)
  {
    assert_eq(num_pcls_sent,0);
    return;
  }

  // velocities are not capped, so all-reduce is needed
  // to check if all particles have been communicated
  //
  long long total_num_pcls_sent = mpi_global_sum(num_pcls_sent);
  if(print_pcl_comm_counts())
  {
    dprintf("spec %d pcls sent: %d, %d",
      is, num_pcls_sent, total_num_pcls_sent);
  }

  if(!total_num_pcls_sent)
    return;

  // apply boundary conditions to incoming particles
  // until they are in the domain to ensure cap on number of
  // communications
  flush_send(); // flush sending of particles
  num_pcls_sent = handle_received_particles(PclCommMode::do_apply_BCs_globally);
  total_num_pcls_sent = mpi_global_sum(num_pcls_sent);
    
  // the maximum number of neighbor communications that would
  // be needed to put a particle in the correct mesh cell
  int comm_max_times = vct->getXLEN()+vct->getYLEN()+vct->getZLEN();
  if(!do_apply_periodic_BC_global) comm_max_times*=2;
  comm_max_times-=3;
  int comm_count=0;
  while(total_num_pcls_sent)
  {
    if(comm_count>=(comm_max_times))
    {
      dprintf("particles still uncommunicated:");
      flush_send();
      num_pcls_sent = handle_received_particles(PclCommMode::print_sent_pcls);
      eprintf("failed to finish up particle communication"
        " within %d communications", comm_max_times);
    }
  
    // flush sending of particles
    flush_send();
    num_pcls_sent = handle_received_particles();
  
    total_num_pcls_sent = mpi_global_sum(num_pcls_sent);
    if(print_pcl_comm_counts())
      dprint(total_num_pcls_sent);
    comm_count++;
  }
}

// exchange particles with neighboring processors
//
// sent particles are deleted from _pcls.
// holes are filled with particles from end.
// then received particles are appended to end.
//
//void Particles3Dcomm::communicate_particles()
//{
//  timeTasks_set_communicating(); // communicating until end of scope
//
//  separate_and_send_particles();
//
//  receive_particles_until_done();
//}

/** return the Kinetic energy */
double Particles3Dcomm::getKe()const
{
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register int i = 0; i < _pcls.size(); i++)
  {
    const SpeciesParticle& pcl = _pcls[i];
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
double Particles3Dcomm::getP()const
{
  double localP = 0.0;
  double totalP = 0.0;
  for (register int i = 0; i < _pcls.size(); i++)
  {
    const SpeciesParticle& pcl = _pcls[i];
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
double Particles3Dcomm::getMaxVelocity()const
{
  double localVel = 0.0;
  double maxVel = 0.0;
  for (int i = 0; i < _pcls.size(); i++)
  {
    const SpeciesParticle& pcl = _pcls[i];
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
// I think that we should be using double rather than long
// long.  This method ignores the weight of the charges.
// This method should be changed to sum the amount of charge
// with each velocity range, which is a double, not an
// integer.  Even in the case that we are counting, double
// precision will work fine.
//
long long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel)const
{
  long long *f = new long long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (int i = 0; i < _pcls.size(); i++) {
    const SpeciesParticle& pcl = _pcls[i];
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
  cout << "Xin = " << beg_pos[0] << "; Xfin = " << end_pos[0] << endl;
  cout << "Yin = " << beg_pos[1] << "; Yfin = " << end_pos[1] << endl;
  cout << "Zin = " << beg_pos[2] << "; Zfin = " << end_pos[2] << endl;
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

void Particles3Dcomm::sort_particles()
{
  timeTasks_begin_task(TimeTasks::MOMENT_PCL_SORTING);
    sort_particles_serial();
  timeTasks_end_task(TimeTasks::MOMENT_PCL_SORTING);
}
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
//      const double xpos = xpos[pidx];
//      const double ypos = ypos[pidx];
//      const double zpos = zpos[pidx];
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

// This can be called from within an omp parallel block
void Particles3Dcomm::copyParticlesToSoA()
{
  const int nop = _pcls.size();
  // create memory for SoA representation
  if(is_output_thread()) dprintf("copying to struct of arrays");
  resize_SoA(nop);
  assert_eq(u.size(),_pcls.size());
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_SOA);
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
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx+=8)
  {
    F64vec8* SoAdata[8] = {
      (F64vec8*) &u[pidx],
      (F64vec8*) &v[pidx],
      (F64vec8*) &w[pidx],
      (F64vec8*) &q[pidx],
      (F64vec8*) &x[pidx],
      (F64vec8*) &y[pidx],
      (F64vec8*) &z[pidx],
      (F64vec8*) &t[pidx]};
    const F64vec8* AoSdata = reinterpret_cast<F64vec8*>(&_pcls[pidx]);
    // this seems to perform slower
    //const F64vec8* AoSdata[8] ={
    //  (F64vec8*) &_pcls[pidx],
    //  (F64vec8*) &_pcls[pidx+1],
    //  (F64vec8*) &_pcls[pidx+2],
    //  (F64vec8*) &_pcls[pidx+3],
    //  (F64vec8*) &_pcls[pidx+4],
    //  (F64vec8*) &_pcls[pidx+5],
    //  (F64vec8*) &_pcls[pidx+6],
    //  (F64vec8*) &_pcls[pidx+7]};
    transpose_8x8_double(SoAdata,AoSdata);
  }
 #endif // __MIC__
  particleType = ParticleType::synched;
}

// This can be called from within an omp parallel block
void Particles3Dcomm::copyParticlesToAoS()
{
  const int nop = u.size();
  if(is_output_thread()) dprintf("copying to array of structs");
  resize_AoS(nop);
  assert_eq(u.size(),_pcls.size());
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_AOS);
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
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx+=8)
  {
    F64vec8* AoSdata = reinterpret_cast<F64vec8*>(&_pcls[pidx]);
    const F64vec8* SoAdata[8] ={
      (F64vec8*) &u[pidx],
      (F64vec8*) &v[pidx],
      (F64vec8*) &w[pidx],
      (F64vec8*) &q[pidx],
      (F64vec8*) &x[pidx],
      (F64vec8*) &y[pidx],
      (F64vec8*) &z[pidx],
      (F64vec8*) &t[pidx]};
    transpose_8x8_double(AoSdata, SoAdata);
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

