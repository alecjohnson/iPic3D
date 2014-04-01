/*******************************************************************************************
  Particles3Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef Part3DCOMM_H
#define Part3DCOMM_H

#include "CollectiveIO.h"
#include "VirtualTopology3D.h"
#include "Grid.h"
#include "Field.h"
#include "Particle.h"
/**
 * 
 * class for particles of the same species with communications methods
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3Dcomm // :public Particles
{
public:
  /** constructor */
  Particles3Dcomm();
  /** destructor */
  ~Particles3Dcomm();
  /** allocate particles */
  void allocate(int species, CollectiveIO * col, VirtualTopology3D * vct, Grid * grid);

  /** calculate the weights given the position of particles */
  //void calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid);
  /** interpolation method GRID->PARTICLE order 1: CIC */
  void interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes */
  int communicate(VirtualTopology3D * ptVCT);
  /** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXleft(double *b_, int np, VirtualTopology3D * vct);
  /** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXright(double *b_, int np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYleft(double *b_, int np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYright(double *b_, int np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZleft(double *b_, int np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZright(double *b_, int np, VirtualTopology3D * vct);
  /** Delete the a particle from a list(array) and pack the list(array) */
  void del_pack(int np, int *nplast);

  /** method to debuild the buffer received */
  int unbuffer(double *b_);

  /** resize the receiving buffer */
  void resize_buffers(int new_buffer_size);
  /** a method to compute how many particles are not in the right domain */
  int isMessagingDone(VirtualTopology3D * ptVCT);
  /** calculate the maximum number exiting from this domain */
  int maxNpExiting();
  /** calculate the weights given the position of particles */
  // void calculateWeights(double*** weight, double xp, double yp, double zp,int ix, int iy, int iz, Grid* grid);

 private:
  void copyParticlesToAoS();
  void copyParticlesToSoA();

 public:
  void convertParticlesToAoS();
  void convertParticlesToSoA();

  /*! sort particles for vectorized push (needs to be parallelized) */
  void sort_particles_serial_SoA_by_xavg(Grid * grid, VirtualTopology3D * vct);
  void sort_particles_serial(Grid * grid, VirtualTopology3D * vct);
  void sort_particles_serial_AoS(Grid * grid, VirtualTopology3D * vct);
  void sort_particles_serial_SoA(Grid * grid, VirtualTopology3D * vct);

  // get accessors for optional arrays
  //
  SpeciesParticle *fetch_pcls(){ return _pcls; }
  SpeciesParticle *fetch_pclstmp(){ return _pclstmp; }
  double * fetch_xavg() { return _xavg; }
  double * fetch_yavg() { return _yavg; }
  double * fetch_zavg() { return _zavg; }
  double * fetch_xtmp() { return _xtmp; }
  double * fetch_ytmp() { return _ytmp; }
  double * fetch_ztmp() { return _ztmp; }
  double * fetch_utmp() { return _utmp; }
  double * fetch_vtmp() { return _vtmp; }
  double * fetch_wtmp() { return _wtmp; }
  double * fetch_qtmp() { return _qtmp; }
  double * fetch_xavgtmp() { return _xavgtmp; }
  double * fetch_yavgtmp() { return _yavgtmp; }
  double * fetch_zavgtmp() { return _zavgtmp; }
  long long *fetch_ParticleIDtmp(){ return _ParticleIDtmp; }

  // inline get accessors
  //
  double get_dx(){return dx;}
  double get_dy(){return dy;}
  double get_dz(){return dz;}
  double get_invdx(){return inv_dx;}
  double get_invdy(){return inv_dy;}
  double get_invdz(){return inv_dz;}
  double get_xstart(){return xstart;}
  double get_ystart(){return ystart;}
  double get_zstart(){return zstart;}
  ParticleType::Type get_particleType()const { return particleType; }
  const SpeciesParticle& get_pcl(int pidx)const{ return _pcls[pidx]; }
  double *getXall()  const { return (x); }
  double *getYall()  const { return (y); }
  double *getZall()  const { return (z); }
  double *getUall()  const { return (u); }
  double *getVall()  const { return (v); }
  double *getWall()  const { return (w); }
  long long *getParticleIDall()  const { return (ParticleID); }
  double *getQall()  const { return (q); }
  // accessors for particle with index indexPart
  double getX(int indexPart)  const { return (x[indexPart]); }
  double getY(int indexPart)  const { return (y[indexPart]); }
  double getZ(int indexPart)  const { return (z[indexPart]); }
  double getU(int indexPart)  const { return (u[indexPart]); }
  double getV(int indexPart)  const { return (v[indexPart]); }
  double getW(int indexPart)  const { return (w[indexPart]); }
  long long getParticleID(int indexPart)  const
    { return (ParticleID[indexPart]); }
  double getQ(int indexPart)  const { return (q[indexPart]); }
  int getNOP()  const { return (nop); }

  // computed get access
  //
  /** return the Kinetic energy */
  double getKe();
  /** return the maximum kinetic energy */
  double getMaxVelocity();
  /** return energy distribution */
  long long *getVelocityDistribution(int nBins, double maxVel);
  /** return the momentum */
  double getP();
  /** Print particles info: positions, velocities */
  void Print(VirtualTopology3D * ptVCT) const;
  /** Print the number of particles of this subdomain */
  void PrintNp(VirtualTopology3D * ptVCT) const;

public:
  // accessors
  int get_ns()const{return ns;}
  int get_numpcls_in_bucket(int cx, int cy, int cz)const
  { return (*numpcls_in_bucket)[cx][cy][cz]; }
  int get_bucket_offset(int cx, int cy, int cz)const
  { return (*bucket_offset)[cx][cy][cz]; }

protected:
  /** number of this species */
  int ns;
  /** maximum number of particles of this species on this domain. used for memory allocation */
  int npmax;
  /** number of particles of this species on this domain */
  int nop;
  /** total number of particles */
  long long np_tot;
  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;
  /** charge to mass ratio */
  double qom;
  /** recon thick */
  double delta;
  /** thermal velocity  - Direction X*/
  double uth;
  /** thermal velocity  - Direction Y*/
  double vth;
  /** thermal velocity  - Direction Z*/
  double wth;
  /** u0 Drift velocity - Direction X */
  double u0;
  /** v0 Drift velocity - Direction Y */
  double v0;
  /** w0 Drift velocity - Direction Z */
  double w0;

  ParticleType::Type particleType;
  // particles data
  //
  // SoA representation
  //
  /** Positions array - X component */
  double *x;
  /** Positions array - Y component */
  double *y;
  /** Positions array - Z component */
  double *z;
  /** Velocities array - X component */
  double *u;
  /** Velocities array - Y component */
  double *v;
  /** Velocities array - Z component */
  double *w;
  /** Charge array */
  double *q;
  /** TrackParticleID */
  bool TrackParticleID;
  /** ParticleID */
  long long *ParticleID;
  //
  // AoS representation
  //
  SpeciesParticle *_pcls;

  // structures for sorting particles
  //
  /** Average position data (used during particle push) **/
  //
  double *_xavg;
  double *_yavg;
  double *_zavg;
  //
  // alternate temporary storage for sorting particles
  //
  long long *_ParticleIDtmp;
  double *_xtmp;
  double *_ytmp;
  double *_ztmp;
  double *_utmp;
  double *_vtmp;
  double *_wtmp;
  double *_qtmp;
  SpeciesParticle *_pclstmp;
  double *_xavgtmp;
  double *_yavgtmp;
  double *_zavgtmp;
  //
  // references for buckets
  //
  array3_int* numpcls_in_bucket;
  array3_int* numpcls_in_bucket_now; // accumulator used during sorting
  //array3_int* bucket_size; // maximum number of particles in bucket
  array3_int* bucket_offset;
  // 
  // bucket totals per thread
  //
  //int num_threads;
  //array3_int* numpcls_in_bucket_thr;
  //arr3_int fetch_numpcls_in_bucket_thr(int i)
  //{
  //  assert_le(0,i);
  //  assert_lt(i,num_threads);
  //  return *(numpcls_in_bucket_thr[i]);
  //};

  /** rank of processor in which particle is created (for ID) */
  int BirthRank[2];
  /** number of variables to be stored in buffer for communication for each particle  */
  int nVar;
  /** time step */
  double dt;
  //
  // Copies of grid data (should just put pointer to Grid in this class)
  //
  /** Simulation domain lengths */
  double xstart, xend, ystart, yend, zstart, zend, invVOL;
  /** Lx = simulation box length - x direction   */
  double Lx;
  /** Ly = simulation box length - y direction   */
  double Ly;
  /** Lz = simulation box length - z direction   */
  double Lz;
  /** grid spacings */
  double dx, dy, dz;
  /** number of grid nodes */
  int nxn, nyn, nzn;
  /** number of grid cells */
  int nxc, nyc, nzc;
  // convenience values from grid
  double inv_dx;
  double inv_dy;
  double inv_dz;
  //
  // Communication variables
  //
  /** buffers for communication */
  /** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
  int buffer_size;
  /** smaller buffer size */
  int buffer_size_small;
  /** buffer with particles going to the right processor - Direction X */
  double *b_X_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_X_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction X */
  double *b_X_LEFT;
  /** pointer to the buffer for resizing */
  double *b_X_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Y */
  double *b_Y_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Y_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Y */
  double *b_Y_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Y_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Z */
  double *b_Z_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Z_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Z */
  double *b_Z_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Z_LEFT_ptr;

  /** number of particles exiting per cycle*/
  int npExitXright;
  /** number of particles exiting to X-LEFT per cycle*/
  int npExitXleft;
  /** number of particles exiting to Y-RIGHT per cycle*/
  int npExitYright;
  /** number of particles exiting to Y-LEFT per cycle*/
  int npExitYleft;
  /** number of particles exiting to Z-RIGHT per cycle*/
  int npExitZright;
  /** number of particles exiting to Z-LEFT per cycle*/
  int npExitZleft;
  /** total number of particles exiting per cycle */
  int npExit;
  /** number of particles not in the right domain   */
  int rightDomain;


  /** bool for communication verbose */
  bool cVERBOSE;
  /** Boundary condition on particles:
          <ul>
          <li>0 = exit</li>
          <li>1 = perfect mirror</li>
          <li>2 = riemission</li>
          <li>3 = periodic condition </li>
          </ul>
          */
  /** Boundary Condition Particles: FaceXright */
  int bcPfaceXright;
  /** Boundary Condition Particles: FaceXleft */
  int bcPfaceXleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceYright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceYleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceZright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceZleft;
  //
  // Other variables
  //
  /** speed of light in vacuum */
  double c;
  /** restart variable for loading particles from restart file */
  int restart;
  /** Number of iteration of the mover*/
  int NiterMover;
  /** velocity of the injection of the particles */
  double Vinj;
  /** removed charge from species */
  double Q_removed;
  /** density of the injection of the particles */
  double Ninj;
};

typedef Particles3Dcomm Particles;

#endif
