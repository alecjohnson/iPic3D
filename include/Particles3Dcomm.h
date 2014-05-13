/*******************************************************************************************
  Particles3Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef Part3DCOMM_H
#define Part3DCOMM_H

//#include "CollectiveIO.h"
// taken from #include "ipicfwd.h" on merge2amaya
class Collective;
typedef Collective CollectiveIO;
#include "Alloc.h"
#include "Particle.h" // for ParticleType
class Grid3DCU;
typedef Grid3DCU Grid;
class EMfields3D;
typedef EMfields3D Field;
class VirtualTopology3D;
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
  Particles3Dcomm(int species, CollectiveIO * col,
    VirtualTopology3D * vct, Grid * grid);
  /** destructor */
  ~Particles3Dcomm();

  /** interpolation method GRID->PARTICLE order 1: CIC */
  // This does not belong in this class and is no longer in use.
  void interpP2G(Field * EMf);

  // communicate particles between processes
  int communicate_particles();

 private:
  void copyParticlesToAoS();
  void copyParticlesToSoA();

 public:
  void convertParticlesToAoS();
  void convertParticlesToSoA();

  /*! sort particles for vectorized push (needs to be parallelized) */
  //void sort_particles_serial_SoA_by_xavg();
  void sort_particles_serial();
  void sort_particles_serial_AoS();
  //void sort_particles_serial_SoA();

  // get accessors for optional arrays
  //
  Larray<SpeciesParticle>& fetch_pcls(){ return _pcls; }
  Larray<SpeciesParticle>& fetch_pclstmp(){ return _pclstmp; }
  Larray<double>& fetch_xavg() { return _xavg; }
  Larray<double>& fetch_yavg() { return _yavg; }
  Larray<double>& fetch_zavg() { return _zavg; }

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
  const SpeciesParticle& get_pcl(int pidx)const{ return fetch_pcls()[pidx]; }
  double *getUall()  const { return &u[0]; }
  double *getVall()  const { return &v[0]; }
  double *getWall()  const { return &w[0]; }
  double *getQall()  const { return &q[0]; }
  double *getXall()  const { return &x[0]; }
  double *getYall()  const { return &y[0]; }
  double *getZall()  const { return &z[0]; }
  //long long *getParticleIDall()  const { return (long long *) q; }
  // accessors for particle with index indexPart
  double getX(int indexPart)  const { return (x[indexPart]); }
  double getY(int indexPart)  const { return (y[indexPart]); }
  double getZ(int indexPart)  const { return (z[indexPart]); }
  double getU(int indexPart)  const { return (u[indexPart]); }
  double getV(int indexPart)  const { return (v[indexPart]); }
  double getW(int indexPart)  const { return (w[indexPart]); }
  double getQ(int indexPart)  const { return (q[indexPart]); }
  //long long& fetch_ParticleID(int indexPart)  const
  //  { return (long long)(q[indexPart]); }
  int getNOP()  const { return (nop); }
  int get_npmax() const {return npmax;}

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
  // velocity components
  Larray<double>& u;
  Larray<double>& v;
  Larray<double>& w;
  // charge
  Larray<double>& q;
  // position
  Larray<double>& x;
  Larray<double>& y;
  Larray<double>& z;
  // subcycle time
  Larray<double>& t;
  // is this class for tracking particles?
  bool TrackParticleID;
  //
  // AoS representation
  //
  Larray<SpeciesParticle>& _pcls;

  // structures for sorting particles
  //
  /** Average position data (used during particle push) **/
  //
  Larray<double>& _xavg;
  Larray<double>& _yavg;
  Larray<double>& _zavg;
  //
  // alternate temporary storage for sorting particles
  //
  Larray<SpeciesParticle>& _pclstmp;
  //
  // references for buckets for serial sort.
  //
  array3_int* numpcls_in_bucket;
  array3_int* numpcls_in_bucket_now; // accumulator used during sorting
  //array3_int* bucket_size; // maximum number of particles in bucket
  array3_int* bucket_offset;

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
  //
  // send buffers
  //
  BlockCommunicator<SpeciesParticle>& sendXleft;
  BlockCommunicator<SpeciesParticle>& sendXrght;
  BlockCommunicator<SpeciesParticle>& sendYleft;
  BlockCommunicator<SpeciesParticle>& sendYrght;
  BlockCommunicator<SpeciesParticle>& sendZleft;
  BlockCommunicator<SpeciesParticle>& sendZrght;
  //
  // recv buffers
  //
  BlockCommunicator<SpeciesParticle>& recvXleft;
  BlockCommunicator<SpeciesParticle>& recvXrght;
  BlockCommunicator<SpeciesParticle>& recvYleft;
  BlockCommunicator<SpeciesParticle>& recvYrght;
  BlockCommunicator<SpeciesParticle>& recvZleft;
  BlockCommunicator<SpeciesParticle>& recvZrght;
  //
  Larray<SpeciesParticle>& bXrght;
  ///** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
  //int buffer_size;
  ///** smaller buffer size */
  //int buffer_size_small;
  ///** buffer with particles going to the right processor - Direction X */
  //double *b_X_RIGHT;
  ///** buffer with particles going to the left processor - Direction X */
  //double *b_X_LEFT;
  ///** buffer with particles going to the right processor - Direction Y */
  //double *b_Y_RIGHT;
  ///** buffer with particles going to the left processor - Direction Y */
  //double *b_Y_LEFT;
  ///** buffer with particles going to the right processor - Direction Z */
  //double *b_Z_RIGHT;
  ///** buffer with particles going to the left processor - Direction Z */
  //double *b_Z_LEFT;

  ///** number of particles exiting per cycle*/
  //int npExitXright;
  ///** number of particles exiting to X-LEFT per cycle*/
  //int npExitXleft;
  ///** number of particles exiting to Y-RIGHT per cycle*/
  //int npExitYright;
  ///** number of particles exiting to Y-LEFT per cycle*/
  //int npExitYleft;
  ///** number of particles exiting to Z-RIGHT per cycle*/
  //int npExitZright;
  ///** number of particles exiting to Z-LEFT per cycle*/
  //int npExitZleft;
  ///** total number of particles exiting per cycle */
  //int npExit;
  ///** number of particles not in the right domain   */
  //int rightDomain;


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
  ///** Boundary Condition Particles: FaceXright */
  //int bcPfaceXright;
  ///** Boundary Condition Particles: FaceXleft */
  //int bcPfaceXleft;
  ///** Boundary Condition Particles: FaceYright */
  //int bcPfaceYright;
  ///** Boundary Condition Particles: FaceYleft */
  //int bcPfaceYleft;
  ///** Boundary Condition Particles: FaceYright */
  //int bcPfaceZright;
  ///** Boundary Condition Particles: FaceYleft */
  //int bcPfaceZleft;
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
