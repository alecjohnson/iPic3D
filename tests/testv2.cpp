#include <omp.h>
#include <stdio.h>
//#include "time.h" // for clock_gettime()
#include <stdint.h> // for uint64_t
#include <stdlib.h> // rand()
#include <sys/time.h>
#include <limits.h> // RAND_MAX
#include <string.h> // memcpy
#include <stdint.h> // for uint64_t etc.
#include <algorithm> // for std::max
#include <cmath> // for std::abs
#include <assert.h>
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"

#define ALIGNMENT 64
#ifdef __INTEL_COMPILER
  #define ALLOC_ALIGNED __attribute__((aligned(ALIGNMENT)));
  #define ASSERT_ALIGNED(X) __assume_aligned(X, ALIGNMENT)
#else
    #define ALLOC_ALIGNED
    #define ASSERT_ALIGNED(X)
#endif

//******** timing ************

#if 0
// measure time using cycle counters
//
inline timespec get_timespec()
{
  // uncomment the desired choice of clock
  //
  // system-wide realtime clock
  //clockid_t clk_id = REALTIME;
  // timer provided by the CPU for each process
  //clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
  // timer provided by the CPU for each thread
  clockid_t clk_id = CLOCK_THREAD_CPUTIME_ID;

  timespec thetime;
  clock_gettime(clk_id, &thetime);
  return thetime;
}

inline timespec diff_timespec(timespec start, timespec end)
{
  timespec diff;
  diff.tv_nsec = end.tv_nsec-start.tv_nsec;
  diff.tv_sec = end.tv_sec-start.tv_sec;
  if(diff.tv_nsec<0) // must we borrow?
  {
    diff.tv_nsec+=1e9;
    diff.tv_sec-=1;
  }
  return diff;
}

inline double get_sec(timespec start_time)
{
  timespec end_time = get_timespec();
  timespec diff_time = diff_timespec(start_time,end_time);
  return diff_time.tv_sec+1.e-9*diff_time.tv_nsec;
}

inline double get_msec(timespec start_time)
{
  return 1.e3*get_sec(start_time);
}
// measure time using clock_gettime
// (more accurate, but requires -lrt when compiling
// and for some reason forces me to compile on
// miclogin rather than knc1 or I can't run the binary.)
//
#define Time timespec
#define get_time() get_timespec()
#define report_time(start) \
 { \
   printf("(in thread %d) %s took %g ms.\n", \
     omp_get_thread_num(), __func__, get_msec(start)); \
 }
#endif


#if 0
// measure wall time
inline double time_msec()
{
  static struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec + tv.tv_usec * 1.e-6)*1.e3;
  // this is another way
  return 1.e3*omp_get_wtime();
}

// measure time using gettimeofday
//
#define Time double
#define report_time(start) \
 { \
   const double end = time_msec(); \
   printf("(in thread %d) %s took %g ms.\n", \
     omp_get_thread_num(), __func__, end - start); \
 }
#define get_time() time_msec()
#endif

// cycle-level timing
//
// read time stamp counter
// could also call __rdtsc intrinsic
// taken from
// http://stackoverflow.com/questions/13772567/get-cpu-cycle-count
//typedef unsigned long uint64_t;
uint64_t rdtsc(){
    unsigned int lo,hi;
    // could first call cpuid to serialize.
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

// measure time using timestamp counter
#define Time uint64_t
#define report_time(start) \
 { \
   uint64_t end = rdtsc(); \
   printf("(in thread %d) %s took %8.6f Mcycles\n", \
     omp_get_thread_num(), __func__, (end - start)*1.e-6); \
 }
 //  printf("(in thread %d) %s took %d cycles\n", \
 //    omp_get_thread_num(), __func__, end - start); \
 //
#define get_time() rdtsc()

//******** tests ************

const int NiterMover=3;
const int D=4;
// const int NPB=2; // number of particles in a block
//#define D (4)
//#define NPB (8)
const int NPB=8;
double field_transpose[6][8] ALLOC_ALIGNED;
double weights[8]ALLOC_ALIGNED;
double Bx;
double By;
double Bz;
double Ex;
double Ey;
double Ez;
double field_components[8][8]ALLOC_ALIGNED;
double fields[8]ALLOC_ALIGNED;

void test_alex()
{
   #pragma simd
   for(int c=0; c<8; c++)
   {
      Bx += weights[c] * field_transpose[0][c];
      By += weights[c] * field_transpose[1][c];
      Bz += weights[c] * field_transpose[2][c];
      Ex += weights[c] * field_transpose[3][c];
      Ey += weights[c] * field_transpose[4][c];
      Ez += weights[c] * field_transpose[5][c];
   }
}

void test_alex2()
{
   for(int i=0; i<8; i++)
   #pragma simd
   for(int c=0; c<8; c++)
   {
      fields[i] += weights[c]*field_transpose[i][c];
   }
}

// unfortunately it seems that I have to make this non-in-lined
// in order to get it to vectorize the way I want.
//
void sample_field(
  double fields[2*D],
  double weights[8],
  double field_components[8][2*D]);

//void sample_field(
//  double fields[2*D],
//  double weights[8],
//  double field_components[8][2*D])__attribute__((noinline));

void test_AoS()
{
   double B[4]ALLOC_ALIGNED;
   double E[4]ALLOC_ALIGNED;
   sample_field(fields,weights,field_components);
   //for(int c=0; c<8; c++)
   //#pragma simd
   //for(int i=0;i<8;i++)
   //{
   //   fields[i] += weights[c]*field_components[c][i];
   //}
   #pragma simd
   for(int i=0;i<4;i++)
   {
     B[i] = fields[i];
     E[i] = fields[4+i];
   }
   for(int i=0;i<4;i++)
   {
     printf("B[%d]=%g\n",i,B[i]);
     printf("E[%d]=%g\n",i,E[i]);
   }
}

// sample from (0,1)
double usample()
{
  // sample from [0,1]
  const double RAND_MAX_inv = 1./RAND_MAX;
  return rand()*RAND_MAX_inv;
  //
  // sample from (0,1)
  const double max_inv = 1./(double(RAND_MAX)+2);
  return (double(rand())+1)*max_inv;
}

// 512 bytes (cache-line-sized)
struct FloatPcl
{
  float x;
  float y;
  float z;
  float t; // time pushed so far
  float u;
  float v;
  float w;
  float q; // charge
  float xavg;
  float yavg;
  float zavg;
  float qom; // charge-to-mass ratio
  float m; // mass
  long long ID;
  int32_t tmp; // filler
};

// 512 bits (cache-line-sized)
struct SpeciesPcl
{
  double x[3];
  double q;
  double u[3];
  long long ID;
 public:
  double get_x(int i){return x[i];}
  double get_u(int i){return u[i];}
  void set_x(int i, double in){x[i]=in;}
  void set_u(int i, double in){u[i]=in;}
  void init_random(int i)
  {
    for(int j=0;j<3;j++)
    {
      x[j] = usample();
      u[j] = usample();
    }
    q=usample();
    ID=i;
  }
};

// 8 particles
class PclBlock
{
  double data[8][8];
 public:
  // convert between AoS and SoA
  // (should be implemented with intrinsics)
  void transpose()
  {
    double temp[8][8]ALLOC_ALIGNED;
    for(int i=0;i<8;i++)
    for(int j=0;j<8;j++)
    {
      temp[i][j] = data[j][i];
    }
    memcpy(&data[0][0],&temp[0][0],sizeof(double)*8*8);
  }
  // assumes AoS
  SpeciesPcl& fetch_pcl(int p){return (SpeciesPcl&)data[p];}
  // assumes SoA
  double* fetch_x(){return data[0];}
  double* fetch_y(){return data[1];}
  double* fetch_z(){return data[2];}
  double* fetch_q(){return data[3];}
  double* fetch_u(){return data[4];}
  double* fetch_v(){return data[5];}
  double* fetch_w(){return data[6];}
  long long* fetch_ID(){return (long long*) data[7];}
};

// 512 bits (cache-line-sized)
class MiPcl
{
  float dx[3];
  float m; // mass
  float xavg[3];
  float qom; // charge-to-mass ratio
  float u[3];
  long long ID;
  int cx[3];
};

double dt;
double dx;
double dy;
double dz;
const int NUMPCLS=NPB*256;
const int NUMBLKS=(NUMPCLS-1)/8+1;
SpeciesPcl pcls[NUMPCLS]ALLOC_ALIGNED;
SpeciesPcl pcls2[NUMPCLS]ALLOC_ALIGNED;
PclBlock pclBlocks[NUMBLKS]ALLOC_ALIGNED;
double x[NUMPCLS]ALLOC_ALIGNED;
double y[NUMPCLS]ALLOC_ALIGNED;
double z[NUMPCLS]ALLOC_ALIGNED;
double u[NUMPCLS]ALLOC_ALIGNED;
double v[NUMPCLS]ALLOC_ALIGNED;
double w[NUMPCLS]ALLOC_ALIGNED;
double q[NUMPCLS]ALLOC_ALIGNED;
long long ID[NUMPCLS]ALLOC_ALIGNED;
// portion of dt remaining for particle
double tpcl[NUMPCLS]ALLOC_ALIGNED;
int destination[NUMPCLS]ALLOC_ALIGNED;

void copy_SoA_to_AoS_particles()
{
  #pragma omp parallel
  {
  const Time start = get_time();
  #pragma omp for
  #pragma simd
  for(int i=0;i<NUMPCLS;i++)
  {
     pcls[i].set_x(0,x[i]);
     pcls[i].set_x(1,y[i]);
     pcls[i].set_x(2,z[i]);
     pcls[i].set_u(0,u[i]);
     pcls[i].set_u(1,v[i]);
     pcls[i].set_u(2,w[i]);
     pcls[i].q = q[i];
     pcls[i].ID = ID[i];
  }
  report_time(start);
  }
}

void copy_AoS_to_SoA_particles()
{
  #pragma omp parallel
  {
  const Time start = get_time();
  #pragma omp for
  #pragma simd
  for(int i=0;i<NUMPCLS;i++)
  {
    x[i] = pcls[i].get_x(0);
    y[i] = pcls[i].get_x(1);
    z[i] = pcls[i].get_x(2);
    u[i] = pcls[i].get_u(0);
    v[i] = pcls[i].get_u(1);
    w[i] = pcls[i].get_u(2);
    q[i] = pcls[i].q;
    ID[i] = pcls[i].ID;
  }
  report_time(start);
  }
}

// Sebastian's revised version
//
void copy_AoS_particles_to_another_array_rev()
{
  uint64_t *p1 = (uint64_t *) pcls;  ASSERT_ALIGNED(p1);
  uint64_t *p2 = (uint64_t *) pcls2; ASSERT_ALIGNED(p2);
  #pragma omp parallel
  {
  //printf("sizeof(SpeciesPcl): %d\n", sizeof(SpeciesPcl));
  const Time start = get_time();
  // With schedule() we make sure that every thread gets aligned
  // consecutive chunk of pcls to copy
  #pragma omp for schedule(static,NUMPCLS*4) nowait
  for(int p=0;p<NUMPCLS*8;p++)
  {
    //ASSERT_ALIGNED(&pcls2[0]);
    //ASSERT_ALIGNED(&pcls[0]);
    //pcls2[p] = pcls[p];
    p2[p] = p1[p];
  }
  report_time(start);
  }
}

void copy_AoS_particles_to_another_array()
{
  #pragma omp parallel
  {
  const Time start = get_time();
  #pragma omp for
  for(int p=0;p<NUMPCLS;p++)
  {
    ASSERT_ALIGNED(&pcls2[0]);
    ASSERT_ALIGNED(&pcls[0]);
    pcls2[p] = pcls[p];
  }
  report_time(start);
  }
}

void copy_AoS_particles_to_another_array_via_memcpy()
{
  const Time start = get_time();
  memcpy(&pcls2[0],&pcls[0],sizeof(pcls[0])*NUMPCLS);
  report_time(start);
}

void copy_AoS_to_AoS_particle_blocks()
{
  #pragma omp parallel
  {
  const Time start = get_time();
  #pragma omp for
  for(int b=0;b<NUMBLKS;b++)
  {
    PclBlock& pclBlock = pclBlocks[b];
    for(int p=0;p<8;p++)
    {
      SpeciesPcl* pcl = &pclBlock.fetch_pcl(p);
      ASSERT_ALIGNED(pcl);
      ASSERT_ALIGNED(&pcls[0]);
      (*pcl) = pcls[b*8+p];
    }
  }
  report_time(start);
  }
}

void copy_AoS_to_AoS_particle_blocks_via_memcpy()
{
  // How to parallelize memcpy?
  const Time start = get_time();
  memcpy(&pclBlocks[0],&pcls[0],sizeof(pcls[0])*NUMPCLS);
  report_time(start);
}

void initialize_data()
{
  for(int c=0;c<8;c++)
  for(int i=0;i<2*D;i++)
    field_components[c][i] = usample();

  // initialize particles
  //
  dt = usample();
  dprint(dt);
  dx = usample();
  dy = usample();
  dz = usample();
  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].init_random(i);
    // ensure that tpcl is no greater than dt
    tpcl[i] = dt*usample();
  }
  // do everything twice in a row to expose cache issues
  //
  // initialize SoA data
  printf("# first pass through data is 5 times slower than second pass:\n");
  copy_AoS_to_SoA_particles();
  copy_AoS_to_SoA_particles();
  // reinitialize AoS data
  printf("# gathering takes similar time as scattering took:\n");
  copy_SoA_to_AoS_particles();
  copy_SoA_to_AoS_particles();
  // initialize particle blocks
  printf("# this copies sequentially and could be done with a single\n");
  printf("# memcpy command, so why is it so much slower?\n");
  copy_AoS_to_AoS_particle_blocks();
  copy_AoS_to_AoS_particle_blocks();
  // same thing but with memcpy
  printf("# this call to memcpy uses only one thread...\n");
  copy_AoS_to_AoS_particle_blocks_via_memcpy();
  copy_AoS_to_AoS_particle_blocks_via_memcpy();
  // initialize pcls2
  printf("# this likewise copies sequentially:\n");
  copy_AoS_particles_to_another_array();
  copy_AoS_particles_to_another_array();
  // same thing but with memcpy
  printf("# and in this case memcpy seems to require about the same time:\n");
  copy_AoS_particles_to_another_array_via_memcpy();
  copy_AoS_particles_to_another_array_via_memcpy();
}

// test pushing all particles in a mesh cell without trying to vectorize
//
void test_push_pcls_in_cell()
{
  double dto2 = usample();
  double qdto2mc = usample();
  double cellstart[D];
  double dx_inv[D];
  for(int i=0;i<3;i++)
  {
    cellstart[i]=usample();
    dx_inv[i]=usample();
  }

  // iterate over all particles in this mesh cell
  //
  const Time start = get_time();
  for(int pi=0;pi<NUMPCLS;pi+=1)
  {
    SpeciesPcl* pcl = &pcls[pi];
    // because SpeciesPcl fits in cache line:
    ASSERT_ALIGNED(pcl);
    double xorig[D]ALLOC_ALIGNED;
    double xavg[D]ALLOC_ALIGNED;
    double uorig[D]ALLOC_ALIGNED;

    // gather position and velocity data from particle block
    //
    for(int i=0;i<D;i++)
    {
      xorig[i] = pcl->get_x(i);
      uorig[i] = pcl->get_u(i);
    }
    //#pragma simd collapse(2)
    for(int i=0;i<D;i++)
    {
      xavg[i] = xorig[i];
    }

    // sample field for this block of particles
    //
    double B[D]ALLOC_ALIGNED;
    double E[D]ALLOC_ALIGNED;
    {
      double fields[D*2]ALLOC_ALIGNED;
      double weights[8]ALLOC_ALIGNED;
      {
        double w[2][D]ALLOC_ALIGNED;
        #pragma simd
        for(int i=0;i<D;i++)
        {
          w[1][i] = dx_inv[i]*(xavg[i]-cellstart[i]);
        }
        #pragma simd
        for(int i=0;i<D;i++)
        {
          w[0][i] = 1.-w[1][i];
        }
        // This can be done in two vectorized
        // multiplications with one swizzle
        // and two shuffles, but can the compiler see that?
        weights[0] = w[0][0]*w[0][1]*w[0][2]; // weight000
        weights[1] = w[0][0]*w[0][1]*w[1][2]; // weight001
        weights[2] = w[0][0]*w[1][1]*w[0][2]; // weight010
        weights[3] = w[0][0]*w[1][1]*w[1][2]; // weight011
        weights[4] = w[1][0]*w[0][1]*w[0][2]; // weight100
        weights[5] = w[1][0]*w[0][1]*w[1][2]; // weight101
        weights[6] = w[1][0]*w[1][1]*w[0][2]; // weight110
        weights[7] = w[1][0]*w[1][1]*w[1][2]; // weight111
      }
      sample_field(fields,weights,field_components);
      #pragma simd
      // scatter field data for vectorized push
      for(int i=0;i<D;i++)
      {
        B[i] = fields[i];
        E[i] = fields[D+i];
      }
    }

    // use sampled field to push particle block
    //
    double uavg[D]ALLOC_ALIGNED;
    {
      double Om[D]ALLOC_ALIGNED;
      double denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      {
        double omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1.0f/float(omsq_p1);
      }
      double ut[D]ALLOC_ALIGNED;
      double udotOm;
      // solve the position equation
      for(int i=0;i<D;i++)
      {
        ut[i] = uorig[i] + qdto2mc*E[i];
      }
      {
        udotOm = ut[0]*Om[0] + ut[1]*Om[1] + ut[2]*Om[2];
      }
      // solve the velocity equation 
      //
      // #pragma simd -- how do I tell it to recognize the swizzle?
      {
        uavg[0] = (ut[0] + (ut[1] * Om[2] - ut[2] * Om[1] + udotOm * Om[0])) * denom;
        uavg[1] = (ut[1] + (ut[2] * Om[0] - ut[0] * Om[2] + udotOm * Om[1])) * denom;
        uavg[2] = (ut[2] + (ut[0] * Om[1] - ut[1] * Om[0] + udotOm * Om[2])) * denom;
      }
      // update average position
      //#pragma simd collapse(2)
      for(int i=0;i<D;i++)
      {
        xavg[i] = xorig[i] + uavg[i] * dto2;
      }
    }

    // update position of particle (assuming this is last iteration)
    {
      for(int i=0;i<D;i++)
      {
        //pcl->set_x(i, xorig[i] + uavg[i]*dt);
        pcl->set_x(i, 2*xavg[i] - xorig[i]);
        pcl->set_u(i, 2*uavg[i] - uorig[i]);
      }
    }
  }
  report_time(start);
}

void test_push_pcls_in_cell_SoA_vectorized()
{
  const int D=3;
  double dto2 = usample();
  double qdto2mc = usample();
  double cellstart[D];
  double dx_inv[D];
  for(int i=0;i<3;i++)
  {
    cellstart[i]=usample();
    dx_inv[i]=usample();
  }
  #pragma omp parallel num_threads(2)
  {

  // iterate over all particles in this mesh cell
  //
  const Time start = get_time();
  #pragma omp for
  #pragma simd
  for(int pi=0;pi<NUMPCLS;pi+=1)
  {
    double xorig[D];
    double xavg[D];
    double uorig[D];

    // get data from particle
    {
      // copy the particle
      xorig[0] = x[pi];
      xorig[1] = y[pi];
      xorig[2] = z[pi];
      uorig[0] = u[pi];
      uorig[1] = v[pi];
      uorig[2] = w[pi];
    }
    for(int i=0;i<D;i++)
    {
      xavg[i] = xorig[i];
    }

    // sample field for this block of particles
    //
    double B[D];
    double E[D];
    {
      double fields[D*2];
      double weights[8];
      {
        double w[2][D];
        for(int i=0;i<D;i++)
        {
          w[1][i] = dx_inv[i]*(xavg[i]-cellstart[i]);
        }
        for(int i=0;i<D;i++)
        {
          w[0][i] = 1.-w[1][i];
        }
        weights[0] = w[0][0]*w[0][1]*w[0][2]; // weight000
        weights[1] = w[0][0]*w[0][1]*w[1][2]; // weight001
        weights[2] = w[0][0]*w[1][1]*w[0][2]; // weight010
        weights[3] = w[0][0]*w[1][1]*w[1][2]; // weight011
        weights[4] = w[1][0]*w[0][1]*w[0][2]; // weight100
        weights[5] = w[1][0]*w[0][1]*w[1][2]; // weight101
        weights[6] = w[1][0]*w[1][1]*w[0][2]; // weight110
        weights[7] = w[1][0]*w[1][1]*w[1][2]; // weight111
      }
      //#pragma unroll
      for(int c=0; c<8; c++)
      //#pragma unroll
      for(int i=0;i<D*2;i++)
      {
        fields[i] += weights[c]*field_components[c][i];
      }
      // scatter field data for vectorized push
      for(int i=0;i<D;i++)
      {
        B[i] = fields[i];
        E[i] = fields[D+i];
      }
    }

    // use sampled field to push particle block
    //
    double uavg[D];
    {
      double Om[D];
      double denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      //const double omsq_p1 = 1.0+(Om[0] * Om[0] + Om[1] * Om[1] + Om[2] * Om[2]);
      //const double denom = 1.0f / float(omsq_p1);
      {
        double omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1.0f/float(omsq_p1);
        //denom = 1.0/omsq_p1;
      }
      double ut[D];
      double udotOm;
      // solve the position equation
      for(int i=0;i<D;i++)
      {
        ut[i] = uorig[i] + qdto2mc*E[i];
      }
      {
        udotOm = ut[0]*Om[0] + ut[1]*Om[1] + ut[2]*Om[2];
      }
      // solve the velocity equation 
      //
      // #pragma simd -- how do I tell it to recognize the swizzle?
      {
        uavg[0] = (ut[0] + (ut[1] * Om[2] - ut[2] * Om[1] + udotOm * Om[0])) * denom;
        uavg[1] = (ut[1] + (ut[2] * Om[0] - ut[0] * Om[2] + udotOm * Om[1])) * denom;
        uavg[2] = (ut[2] + (ut[0] * Om[1] - ut[1] * Om[0] + udotOm * Om[2])) * denom;
      }
      // update average position
      for(int i=0;i<D;i++)
      {
        xavg[i] = xorig[i] + uavg[i] * dto2;
      }
    }

    // update position of particle (assuming this is last iteration)
    {
      x[pi] = 2*xavg[0] - xorig[0];
      y[pi] = 2*xavg[1] - xorig[1];
      z[pi] = 2*xavg[2] - xorig[2];
      u[pi] = 2*uavg[0] - uorig[0];
      v[pi] = 2*uavg[1] - uorig[1];
      w[pi] = 2*uavg[2] - uorig[2];
    }
  }
  report_time(start);
  }
}

// This should require the same execution time as the previous
// method but in fact takes less time, especially if denom is
// calculated with double precision.
//
void move_bucket_old()
{
    const double dto2 = usample();
    const double qdto2mc = usample();
    const double cellstartx=usample();
    const double cellstarty=usample();
    const double cellstartz=usample();
    const double dx_inv=usample();
    const double dy_inv=usample();
    const double dz_inv=usample();
    #pragma omp parallel 
    {
        const Time start = get_time();
        #pragma omp for
        #pragma simd
        for(int pidx = 0; pidx < NUMPCLS; pidx++)
        {
            // copy the particle
            const double xorig = x[pidx];
            const double yorig = y[pidx];
            const double zorig = z[pidx];
            const double uorig = u[pidx];
            const double vorig = v[pidx];
            const double worig = w[pidx];

            // initialize xavg to xorig
            double xavg = x[pidx];
            double yavg = y[pidx];
            double zavg = z[pidx];

            // compute weights for field components
            //
            double weights[8];
            // fraction of the distance from the left of the cell
            const double w0x = dx_inv*(xavg - cellstartx);
            const double w0y = dy_inv*(yavg - cellstarty);
            const double w0z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the right
            const double w1x = 1.-w0x;
            const double w1y = 1.-w0y;
            const double w1z = 1.-w0z;
            const double weight00 = w0x*w0y;
            const double weight01 = w0x*w1y;
            const double weight10 = w1x*w0y;
            const double weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            double Exl = 0.0;
            double Eyl = 0.0;
            double Ezl = 0.0;
            double Bxl = 0.0;
            double Byl = 0.0;
            double Bzl = 0.0;

            // would expanding this out help to vectorize?
            for(int c=0; c<8; c++)
            {
                Bxl += weights[c] * field_components[c][0];
                Byl += weights[c] * field_components[c][1];
                Bzl += weights[c] * field_components[c][2];
                Exl += weights[c] * field_components[c][3];
                Eyl += weights[c] * field_components[c][4];
                Ezl += weights[c] * field_components[c][5];
            }

            const double Omx = qdto2mc*Bxl;
            const double Omy = qdto2mc*Byl;
            const double Omz = qdto2mc*Bzl;

            // end interpolation
            const double omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const double denom = 1.0f / float(omsq_p1);
            // solve the position equation
            const double ut = uorig + qdto2mc * Exl;
            const double vt = vorig + qdto2mc * Eyl;
            const double wt = worig + qdto2mc * Ezl;
            //const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const double udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const double uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const double vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const double wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
            // update average position
            xavg = xorig + uavg * dto2;
            yavg = yorig + vavg * dto2;
            zavg = zorig + wavg * dto2;

            // update particle (assuming this is last iteration)
            {
                //x[pidx] = xorig + uavg * dt;
                //y[pidx] = yorig + vavg * dt;
                //z[pidx] = zorig + wavg * dt;
                x[pidx] = 2.0 * xavg - xorig;
                y[pidx] = 2.0 * yavg - yorig;
                z[pidx] = 2.0 * zavg - zorig;
                u[pidx] = 2.0 * uavg - uorig;
                v[pidx] = 2.0 * vavg - vorig;
                w[pidx] = 2.0 * wavg - worig;
            }
        }
        report_time(start);
    }
}

// assumes dXfull, dYfull, dZfull are nonnegative
//
// returns desingularized version of
//   min(dXwant/dXfull, dYwant/dYfull, dZwant/dZfull)
// calculated with a single reciprocal.
//
static inline float get_truncation_ratio(
  //int&direction,
  float dXwant, float dXfull,
  float dYwant, float dYfull,
  float dZwant, float dZfull)
{
  const float motion_freedom = 1.e-2;
  // This modification of the input modifies the truncation ratio
  // and is designed to ensure the following properties:
  // * particle is stopped before going
  //   motion_freedom beyond the boundary
  // * if want>full then (modified) ratio > 1
  //   (particles that should not be stopped aren't)
  // * if want<full then modified_ratio > ratio
  //   (so particle is always allowed to leave)
  //
  // If we are willing to accept that particles never stop
  // exactly at cell boundaries, we could replace "max" with
  // addition here.
  //
  // make sure that distance to wall is strictly positive
  //
  dXwant = std::max(motion_freedom,dXwant);
  dYwant = std::max(motion_freedom,dYwant);
  dZwant = std::max(motion_freedom,dZwant);
  //
  // avoid division singularities
  //
  dXfull = std::max(motion_freedom,dXfull);
  dYfull = std::max(motion_freedom,dYfull);
  dZfull = std::max(motion_freedom,dZfull);

  const float denominator = dXfull*dYfull*dZfull;
  const float dXprod = dXwant*dYfull*dZfull;
  const float dYprod = dYwant*dXfull*dZfull;
  const float dZprod = dZwant*dXfull*dYfull;
  float numerator; // = denominator;
  if(dXprod<dYprod)
  {
    if(dXprod<dZprod)
    {
      numerator = dXprod;
      //direction = 1;
    }
    else
    {
      numerator = dZprod;
      //direction = 4;
    }
  }
  else
  {
    if(dYprod<dZprod)
    {
      numerator = dYprod;
      //direction = 2;
    }
    else
    {
      numerator = dZprod;
      //direction = 4;
    }
  }
  return numerator/denominator;
}

// in which direction are we most outside?
// +/-1: +/-X
// +/-2: +/-Y
// +/-4: +/-Z
static inline int get_direction(float Xpos, float Ypos, float Zpos)
{
  int direction;
  const float aX = std::abs(Xpos);
  const float aY = std::abs(Ypos);
  const float aZ = std::abs(Zpos);
  //
  // This way requires 8 comparisons and 7 masked assignments
  //
  //if(aX > aY) direction = (aX>aZ) ? 1 : 4;
  //else direction = (aY>aZ) ? 2 : 4;
  //const float low = max(xpos,ypos,zpos);
  //const float hgh = min(xpos,ypos,zpos);
  //if(hgh < -low) direction = -direction;
  //
  // This way requires 7 comparisons and 8 masked assignments
  //
  if(aX>aY)
  {
    if (aX>aZ)
    {
      if(Xpos > 0)
        direction = 1;
      else
        direction = -1;
    }
    else
    {
      if(Zpos > 0)
        direction = 4;
      else
        direction = -4;
    }
  }
  else
  {
    if(aY>aZ)
    {
      if(Ypos > 0)
        direction = 2;
      else
        direction = -2;
    }
    else
    {
      if(Zpos > 0)
        direction = 4;
      else
        direction = -4;
    }
  }
  return direction;
}

// try to take a time step that stops particle at a face
//
// alternative approach: allow to move at most one mesh cell.
// requires computing norm of motion.  Use infinity norm.
//
void push_pcls_in_cell_SoA_stopping_at_face()
{
    // time step resolution (analogous to FLT_MIN,
    // defines a limit on precision)
    const float dt_min = 1e-6*dt;
    // shortest allowed subcycle time step
    // (assumed to be no less than dt_min)
    const float min_dt = 1e-5*dt;
    //const double dto2 = usample();
    //const double qdto2mc = usample();
    const double qo2mc = usample();
    const double dx_over_two = dx/2;
    const double dy_over_two = dy/2;
    const double dz_over_two = dz/2;
    const double two_over_dx = 2/dx;
    const double two_over_dy = 2/dy;
    const double two_over_dz = 2/dz;
    const double xmiddle = usample(); // position of middle of cell
    const double ymiddle = usample(); // position of middle of cell
    const double zmiddle = usample(); // position of middle of cell
    #pragma omp parallel 
    {
        const Time start = get_time();
        #pragma omp for
        #pragma simd
        for(int pidx = 0; pidx < NUMPCLS; pidx++)
        {
          // copy the particle
          //
          // x is physical position
          // X is position is in canonical coordinates (-1 <=~ X <=~ 1)
          //
          //const float Xorig = X[pidx];
          //const float Yorig = Y[pidx];
          //const float Zorig = Z[pidx];
          const float Xorig = (x[pidx]-xmiddle)*two_over_dx;
          const float Yorig = (y[pidx]-ymiddle)*two_over_dy;
          const float Zorig = (z[pidx]-zmiddle)*two_over_dz;
          // u is physical velocity
          const float uorig = u[pidx];
          const float vorig = v[pidx];
          const float worig = w[pidx];
          const float tpcl_ = tpcl[pidx];
          //assert_le(tpcl_, dt);

          // The computed time step needs to be double
          // precision up to the point where the calculation is
          // unique for every particle.
          //
          // compute time remaining for particle until
          // next synchonization point.  Note that if tpcl_==0
          // then there is no loss of precision at this point.
          double dtpcl = dt-tpcl_;
          // initialize subcycle time to be remaining time
          // (used in first iteration of iterative solver)
          //
          double dtcycle = dtpcl;

          // initialize xavg to xorig
          float Xavg = Xorig;
          float Yavg = Yorig;
          float Zavg = Zorig;

          // purpose of iterative solver is to find
          // dtcycle, Xavg.., and uavg...
          float uavg, vavg, wavg;

          // this is the part that must vectorize
          // #pragma omp simd
          for(int niter=0;niter<NiterMover;niter++)
          {
            // compute weights for field components
            //
            float weights[8];
            // fraction of the distance from the left of the cell
            const float w0x = 0.5*Xavg + 0.5;
            const float w0y = 0.5*Yavg + 0.5;
            const float w0z = 0.5*Zavg + 0.5;
            // fraction of distance from the right
            const double w1x = 1.-w0x;
            const double w1y = 1.-w0y;
            const double w1z = 1.-w0z;
            const double weight00 = w0x*w0y;
            const double weight01 = w0x*w1y;
            const double weight10 = w1x*w0y;
            const double weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            float Exl = 0.0;
            float Eyl = 0.0;
            float Ezl = 0.0;
            float Bxl = 0.0;
            float Byl = 0.0;
            float Bzl = 0.0;

            // field_components is double precision; loss of
            // precision at this point is expected to mitigated
            // by the large number of particles.  When we sum
            // moments we will go back to double precision.
            //
            for(int c=0; c<8; c++)
            {
                Bxl += weights[c] * field_components[c][0];
                Byl += weights[c] * field_components[c][1];
                Bzl += weights[c] * field_components[c][2];
                Exl += weights[c] * field_components[c][3];
                Eyl += weights[c] * field_components[c][4];
                Ezl += weights[c] * field_components[c][5];
            }

            const double qdto2mc = qo2mc*dtcycle;
            const float Omx = qdto2mc*Bxl;
            const float Omy = qdto2mc*Byl;
            const float Omz = qdto2mc*Bzl;

            // end interpolation
            const float omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const float denom = 1.0f / omsq_p1;
            // solve the position equation
            const float ut = uorig + qdto2mc * Exl;
            const float vt = vorig + qdto2mc * Eyl;
            const float wt = worig + qdto2mc * Ezl;
            //const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const float udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const float uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const float vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const float wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

            // stop the particle at the cell boundary
            //
            // compute the displacement assuming the particle is not stopped
            //
            const float dxpcl = dtcycle*uavg;
            const float dypcl = dtcycle*vavg;
            const float dzpcl = dtcycle*wavg;
            const float dXpcl = dxpcl*two_over_dx;
            const float dYpcl = dypcl*two_over_dy;
            const float dZpcl = dzpcl*two_over_dz;
            const float Xnew = Xorig + dXpcl;
            const float Ynew = Yorig + dYpcl;
            const float Znew = Zorig + dZpcl;
            //
            // compute the factor by which the motion
            // (time step) must be multiplied for the particle
            // to stop at the cell boundary.
            //
            // 1. compute the distance moved.
            //
            const float dXmag = std::abs(dXpcl);
            const float dYmag = std::abs(dYpcl);
            const float dZmag = std::abs(dZpcl);
            //
            // 2. compute the distance to the wall
            //    (if moving away then allow no motion)
            //
            float dXwall, dYwall, dZwall;
            const float hghXwall=1., lowXwall=-1.;
            const float hghYwall=1., lowYwall=-1.;
            const float hghZwall=1., lowZwall=-1.;
            // calculate (signed) distance to wall in direction of motion
            if(dXpcl > 0)
              dXwall = hghXwall - Xorig;
            else
              dXwall = Xorig - lowXwall;
            //
            if(dYpcl > 0)
              dYwall = hghYwall - Yorig;
            else
              dYwall = Yorig - lowYwall;
            //
            if(dZpcl > 0)
              dZwall = hghZwall - Zorig;
            else
              dZwall = Zorig - lowZwall;
            //
            // The following alternative avoids comparisons, but
            // unfortunately a singularity arises if the particle
            // begins just outside the cell and is heading toward
            // the cell; dealing with this singularity requires
            // comparisons...
            //
            // distance from orig pos to wall
            // (unless motion would not even bring particle to
            // the middle of the cell, in which case this is
            // the distance to the reflection of the other wall
            // across the final position, which would at least
            // double the distance and therefore should not cause
            // a problem in the final result; a subsequent iteration
            // should then be able to take the particle to the wall...).
            //dXwall = dXmag + 1 - abs(Xnew);
            //dXwall = max(abs(Xorig)-1, dXwall);
            //
            // 3. compute the ratio to truncate motion
            //
            const float ratio = get_truncation_ratio(
              dXwall, dXmag,
              dYwall, dYmag,
              dZwall, dZmag);
            //
            // 4. truncate or adjust the time step and motion accordingly
            //
            const float dtproposed = dtcycle*ratio;
            // enforce a minimum dtcycle
            dtcycle = std::max(min_dt, dtproposed);
            // but cap final time at synchronization time
            // (at which point double precision is here restored)
            dtcycle = std::min(dtpcl, dtcycle);

            // apply the corrected time step
            Xavg = Xorig + 0.5*dtcycle*uavg;
            Yavg = Yorig + 0.5*dtcycle*vavg;
            Zavg = Zorig + 0.5*dtcycle*wavg;
          }

          // update particle after the last iteration
          const float Xend = 2.0 * Xavg - Xorig;
          const float Yend = 2.0 * Yavg - Yorig;
          const float Zend = 2.0 * Zavg - Zorig;
          {
            //X[pidx] = Xnew;
            //Y[pidx] = Ynew;
            //Z[pidx] = Znew;
            x[pidx] = Xend*dx_over_two+xmiddle;
            y[pidx] = Yend*dy_over_two+ymiddle;
            z[pidx] = Zend*dz_over_two+zmiddle;
            u[pidx] = 2.0 * uavg - uorig;
            v[pidx] = 2.0 * vavg - vorig;
            w[pidx] = 2.0 * wavg - worig;
            tpcl[pidx] += dtcycle;
          }

          // if the particle is done being moved then it
          // can be presumed inside the box
          if(tpcl[pidx] > (dt-dt_min))
          {
            destination[pidx] = 0;
          }
          else
          // particle is not done, so move to neighbor
          {
            // for non-canonical coordinates
            //destination[pidx] = get_direction(
            //  (x[pidx]-xmiddle)*two_over_dx,
            //  (y[pidy]-ymiddle)*two_over_dy,
            //  (z[pidz]-zmiddle)*two_over_dx);

            // for canonical coordinates
            destination[pidx] = get_direction(Xend, Yend, Zend);
          }
        }
        report_time(start);
    }
}

// move particles in 8-particle transposable blocks
void move_SoA_blocks()
{
    const double dto2 = usample();
    const double qdto2mc = usample();
    const double cellstartx=usample();
    const double cellstarty=usample();
    const double cellstartz=usample();
    const double dx_inv=usample();
    const double dy_inv=usample();
    const double dz_inv=usample();
    #pragma omp parallel 
    {
      const Time start = get_time();
      #pragma omp for
      for(int bidx = 0; bidx < NUMBLKS; bidx++)
      {
        PclBlock& pclBlock = pclBlocks[bidx];
        pclBlock.transpose();
        double* x = pclBlock.fetch_x();
        double* y = pclBlock.fetch_y();
        double* z = pclBlock.fetch_z();
        double* u = pclBlock.fetch_u();
        double* v = pclBlock.fetch_v();
        double* w = pclBlock.fetch_w();
        ASSERT_ALIGNED(x);
        ASSERT_ALIGNED(y);
        ASSERT_ALIGNED(z);
        ASSERT_ALIGNED(u);
        ASSERT_ALIGNED(v);
        ASSERT_ALIGNED(w);
        #pragma omp simd
        for(int pidx = 0; pidx < 8; pidx++)
        {
            // copy the particle
            const double xorig = x[pidx];
            const double yorig = y[pidx];
            const double zorig = z[pidx];
            const double uorig = u[pidx];
            const double vorig = v[pidx];
            const double worig = w[pidx];

            // initialize xavg to xorig
            double xavg = x[pidx];
            double yavg = y[pidx];
            double zavg = z[pidx];

            // compute weights for field components
            //
            double weights[8];
            // fraction of the distance from the left of the cell
            const double w0x = dx_inv*(xavg - cellstartx);
            const double w0y = dy_inv*(yavg - cellstarty);
            const double w0z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the right
            const double w1x = 1-w0x;
            const double w1y = 1-w0y;
            const double w1z = 1-w0z;
            const double weight00 = w0x*w0y;
            const double weight01 = w0x*w1y;
            const double weight10 = w1x*w0y;
            const double weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            double Exl = 0.0;
            double Eyl = 0.0;
            double Ezl = 0.0;
            double Bxl = 0.0;
            double Byl = 0.0;
            double Bzl = 0.0;

            // would expanding this out help to vectorize?
            for(int c=0; c<8; c++)
            {
                Bxl += weights[c] * field_components[c][0];
                Byl += weights[c] * field_components[c][1];
                Bzl += weights[c] * field_components[c][2];
                Exl += weights[c] * field_components[c][3];
                Eyl += weights[c] * field_components[c][4];
                Ezl += weights[c] * field_components[c][5];
            }

            const double Omx = qdto2mc*Bxl;
            const double Omy = qdto2mc*Byl;
            const double Omz = qdto2mc*Bzl;

            // end interpolation
            const double omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const double denom = 1.0f / float(omsq_p1);
            // solve the position equation
            const double ut = uorig + qdto2mc * Exl;
            const double vt = vorig + qdto2mc * Eyl;
            const double wt = worig + qdto2mc * Ezl;
            //const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const double udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const double uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const double vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const double wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
            // update average position
            xavg = xorig + uavg * dto2;
            yavg = yorig + vavg * dto2;
            zavg = zorig + wavg * dto2;

            // update particle (assuming this is last iteration)
            {
                //x[pidx] = xorig + uavg * dt;
                //y[pidx] = yorig + vavg * dt;
                //z[pidx] = zorig + wavg * dt;
                x[pidx] = 2.0 * xavg - xorig;
                y[pidx] = 2.0 * yavg - yorig;
                z[pidx] = 2.0 * zavg - zorig;
                u[pidx] = 2.0 * uavg - uorig;
                v[pidx] = 2.0 * vavg - vorig;
                w[pidx] = 2.0 * wavg - worig;
            }
        }
        pclBlock.transpose();
      }
      report_time(start);
    }
}

void test_push_pcls_in_cell_AoS_scatter_gather_vectorization()
{
  double dto2 = usample();
  double qdto2mc = usample();
  double cellstart[D];
  double dx_inv[D];
  for(int i=0;i<3;i++)
  {
    cellstart[i]=usample();
    dx_inv[i]=usample();
  }
  #pragma omp parallel num_threads(2)
  {

  // iterate over all particles in this mesh cell
  //
  const Time start = get_time();
  //simd accelerates execution by a factor of 3
  #pragma omp for
  #pragma simd
  for(int pi=0;pi<NUMPCLS;pi+=1)
  {
    SpeciesPcl* pcl = &pcls[pi];
    // because SpeciesPcl fits in cache line:
    ASSERT_ALIGNED(pcl);
    double xorig[D]ALLOC_ALIGNED;
    double xavg[D]ALLOC_ALIGNED;
    double uorig[D]ALLOC_ALIGNED;

    // gather position and velocity data from particle block
    //
    for(int i=0;i<D;i++)
    {
      xorig[i] = pcl->get_x(i);
      uorig[i] = pcl->get_u(i);
    }
    for(int i=0;i<D;i++)
    {
      xavg[i] = xorig[i];
    }

    // sample field for this block of particles
    //
    double B[D]ALLOC_ALIGNED;
    double E[D]ALLOC_ALIGNED;
    {
      double fields[D*2]ALLOC_ALIGNED;
      double weights[8]ALLOC_ALIGNED;
      if(true)
      {
        const double w1x = dx_inv[0]*(xavg[0] - cellstart[0]);
        const double w1y = dx_inv[1]*(xavg[1] - cellstart[1]);
        const double w1z = dx_inv[2]*(xavg[2] - cellstart[2]);
        const double w0x = 1-w1x;
        const double w0y = 1-w1y;
        const double w0z = 1-w1z;
        const double weight00 = w0x*w0y;
        const double weight01 = w0x*w1y;
        const double weight10 = w1x*w0y;
        const double weight11 = w1x*w1y;
        weights[0] = weight00*w0z; // weight000
        weights[1] = weight00*w1z; // weight001
        weights[2] = weight01*w0z; // weight010
        weights[3] = weight01*w1z; // weight011
        weights[4] = weight10*w0z; // weight100
        weights[5] = weight10*w1z; // weight101
        weights[6] = weight11*w0z; // weight110
        weights[7] = weight11*w1z; // weight111
      }
      else
      {
        double w[2][D]ALLOC_ALIGNED;
        for(int i=0;i<D;i++)
        {
          w[1][i] = dx_inv[i]*(xavg[i]-cellstart[i]);
        }
        for(int i=0;i<D;i++)
        {
          w[0][i] = 1.-w[1][i];
        }
        weights[0] = w[0][0]*w[0][1]*w[0][2]; // weight000
        weights[1] = w[0][0]*w[0][1]*w[1][2]; // weight001
        weights[2] = w[0][0]*w[1][1]*w[0][2]; // weight010
        weights[3] = w[0][0]*w[1][1]*w[1][2]; // weight011
        weights[4] = w[1][0]*w[0][1]*w[0][2]; // weight100
        weights[5] = w[1][0]*w[0][1]*w[1][2]; // weight101
        weights[6] = w[1][0]*w[1][1]*w[0][2]; // weight110
        weights[7] = w[1][0]*w[1][1]*w[1][2]; // weight111
      }
      //#pragma unroll
      for(int c=0; c<8; c++)
      //#pragma unroll
      for(int i=0;i<D*2;i++)
      {
        fields[i] += weights[c]*field_components[c][i];
      }
      // scatter field data for vectorized push
      for(int i=0;i<D;i++)
      {
        B[i] = fields[i];
        E[i] = fields[D+i];
      }
    }

    // use sampled field to push particle block
    //
    double uavg[D]ALLOC_ALIGNED;
    {
      double Om[D]ALLOC_ALIGNED;
      double denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      {
        double omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1.0f/float(omsq_p1);
      }
      double ut[D]ALLOC_ALIGNED;
      double udotOm;
      // solve the position equation
      for(int i=0;i<D;i++)
      {
        ut[i] = uorig[i] + qdto2mc*E[i];
      }
      {
        udotOm = ut[0]*Om[0] + ut[1]*Om[1] + ut[2]*Om[2];
      }
      // solve the velocity equation 
      //
      // #pragma simd -- how do I tell it to recognize the swizzle?
      {
        uavg[0] = (ut[0] + (ut[1] * Om[2] - ut[2] * Om[1] + udotOm * Om[0])) * denom;
        uavg[1] = (ut[1] + (ut[2] * Om[0] - ut[0] * Om[2] + udotOm * Om[1])) * denom;
        uavg[2] = (ut[2] + (ut[0] * Om[1] - ut[1] * Om[0] + udotOm * Om[2])) * denom;
      }
      // update average position
      for(int i=0;i<D;i++)
      {
        xavg[i] = xorig[i] + uavg[i] * dto2;
      }
    }

    // update position of particle (assuming this is last iteration)
    {
      for(int i=0;i<D;i++)
      {
        pcl->set_x(i, 2*xavg[i] - xorig[i]);
        pcl->set_u(i, 2*uavg[i] - uorig[i]);
      }
    }
  }
  report_time(start);
  }
}

// test pushing all particles in a mesh cell
void test_push_pcls_in_cell_AoS_localized_vectorization()
{
  double dto2 = usample();
  double qdto2mc = usample();
  double cellstart[D];
  double dx_inv[D];
  for(int i=0;i<3;i++)
  {
    cellstart[i]=usample();
    dx_inv[i]=usample();
  }
  #pragma omp parallel num_threads(2)
  {

  // iterate over all particles in this mesh cell
  //
  const Time start = get_time();
  #pragma omp for
  for(int pi=0;pi<NUMPCLS;pi+=NPB)
  {
    SpeciesPcl* pcl = &pcls[pi];
    // because SpeciesPcl fits in cache line:
    ASSERT_ALIGNED(pcl);
    double xorig[NPB][D]ALLOC_ALIGNED;
    double xavg[NPB][D]ALLOC_ALIGNED;
    double uorig[NPB][D]ALLOC_ALIGNED;

    // gather position and velocity data from particle block
    //
    for(int p=0; p<NPB; p++)
    for(int i=0;i<D;i++)
    {
      xorig[p][i] = pcl[p].get_x(i);
      uorig[p][i] = pcl[p].get_u(i);
    }
    //#pragma simd collapse(2) // this generates 56 instructions
    for(int p=0; p<NPB; p++)
    for(int i=0;i<D;i++)
    {
      xavg[p][i] = xorig[p][i];
    }

    // sample field for this block of particles
    //
    double B[NPB][D]ALLOC_ALIGNED;
    double E[NPB][D]ALLOC_ALIGNED;
    double ws[2][NPB][D]ALLOC_ALIGNED;
    bool vectorized_w = false;
    if(vectorized_w)
    {
      for(int p=0; p<NPB; p++)
      for(int i=0;i<D;i++)
      {
        ws[1][p][i] = dx_inv[i]*(xavg[p][i]-cellstart[i]);
        ws[0][p][i] = 1.-ws[1][p][i];
      }
      double fields[NPB][2*D]ALLOC_ALIGNED;
      for(int p=0; p<NPB; p++)
      {
        double weights[8]ALLOC_ALIGNED;
        {
          // This can be done in two vectorized
          // multiplications with one swizzle
          // and two shuffles, but can the compiler see that?
          weights[0] = ws[0][p][0]*ws[0][p][1]*ws[0][p][2]; // weight000
          weights[1] = ws[0][p][0]*ws[0][p][1]*ws[1][p][2]; // weight001
          weights[2] = ws[0][p][0]*ws[1][p][1]*ws[0][p][2]; // weight010
          weights[3] = ws[0][p][0]*ws[1][p][1]*ws[1][p][2]; // weight011
          weights[4] = ws[1][p][0]*ws[0][p][1]*ws[0][p][2]; // weight100
          weights[5] = ws[1][p][0]*ws[0][p][1]*ws[1][p][2]; // weight101
          weights[6] = ws[1][p][0]*ws[1][p][1]*ws[0][p][2]; // weight110
          weights[7] = ws[1][p][0]*ws[1][p][1]*ws[1][p][2]; // weight111
        }
        sample_field(fields[p],weights,field_components);
        #pragma simd
        // scatter field data for vectorized push
        for(int i=0;i<D;i++)
        {
          B[p][i] = fields[p][i];
          E[p][i] = fields[p][D+i];
        }
      }
    }
    else
    {
      for(int p=0; p<NPB; p++)
      {
        double fields[2*D]ALLOC_ALIGNED;
        double weights[8]ALLOC_ALIGNED;
        {
          double w[2][D]ALLOC_ALIGNED;
          for(int i=0;i<D;i++)
          {
            w[1][i] = dx_inv[i]*(xavg[p][i]-cellstart[i]);
          }
          for(int i=0;i<D;i++)
          {
            w[0][i] = 1.-w[1][i];
          }
          // This can be done in two vectorized
          // multiplications with one swizzle
          // and two shuffles, but can the compiler see that?
          weights[0] = w[0][0]*w[0][1]*w[0][2]; // weight000
          weights[1] = w[0][0]*w[0][1]*w[1][2]; // weight001
          weights[2] = w[0][0]*w[1][1]*w[0][2]; // weight010
          weights[3] = w[0][0]*w[1][1]*w[1][2]; // weight011
          weights[4] = w[1][0]*w[0][1]*w[0][2]; // weight100
          weights[5] = w[1][0]*w[0][1]*w[1][2]; // weight101
          weights[6] = w[1][0]*w[1][1]*w[0][2]; // weight110
          weights[7] = w[1][0]*w[1][1]*w[1][2]; // weight111
        }
        sample_field(fields,weights,field_components);
        // scatter field data for vectorized push
        for(int i=0;i<D;i++)
        {
          B[p][i] = fields[i];
          E[p][i] = fields[D+i];
        }
      }
    }

    // use sampled field to push particle block
    //
    double uavg[NPB][D]ALLOC_ALIGNED;
    {
      double Om[NPB][D]ALLOC_ALIGNED;
      double denom[NPB]ALLOC_ALIGNED;
      for(int p=0;p<NPB;p++)
      for(int i=0;i<D;i++)
      {
        Om[p][i] = qdto2mc*B[p][i];
      }
      for(int p=0;p<NPB;p++)
      {
        double omsq_p1 = 1.
                     + Om[p][0] * Om[p][0]
                     + Om[p][1] * Om[p][1]
                     + Om[p][2] * Om[p][2];
        // This generates 29 instructions
        //denom[p] = 1.0f/float(omsq_p1);
        // This generates 42 instructions
        //denom[p] = 1.0f/omsq_p1;
        // This generates 42 instructions
        denom[p] = 1.0/omsq_p1;
      }
      double ut[NPB][D]ALLOC_ALIGNED;
      double udotOm[NPB]ALLOC_ALIGNED;
      // solve the position equation
      for(int p=0;p<NPB;p++)
      for(int i=0;i<D;i++)
      {
        ut[p][i] = uorig[p][i] + qdto2mc*E[p][i];
      }
      for(int p=0;p<NPB;p++)
      {
        udotOm[p] = ut[p][0]*Om[p][0] + ut[p][1]*Om[p][1] + ut[p][2]*Om[p][2];
      }
      // solve the velocity equation 
      //
      //#pragma simd // how do I tell it to recognize the swizzle?
      for(int p=0;p<NPB;p++)
      {
        uavg[p][0] = (ut[p][0] + (ut[p][1] * Om[p][2] - ut[p][2] * Om[p][1] + udotOm[p] * Om[p][0])) * denom[p];
        uavg[p][1] = (ut[p][1] + (ut[p][2] * Om[p][0] - ut[p][0] * Om[p][2] + udotOm[p] * Om[p][1])) * denom[p];
        uavg[p][2] = (ut[p][2] + (ut[p][0] * Om[p][1] - ut[p][1] * Om[p][0] + udotOm[p] * Om[p][2])) * denom[p];
      }
      // update average position
      for(int p=0;p<NPB;p++)
      for(int i=0;i<D;i++)
      {
        xavg[p][i] = xorig[p][i] + uavg[p][i] * dto2;
      }
    }

    // update position of particle (assuming this is last iteration)
    {
      for(int p=0;p<NPB;p++)
      #pragma simd
      for(int i=0;i<D;i++)
      {
        pcl[p].set_x(i, 2*xavg[p][i] - xorig[p][i]);
        pcl[p].set_u(i, 2*uavg[p][i] - uorig[p][i]);
      }
    }
  }
  report_time(start);
  }
}

// This is declared noinline in order to force vectorization the way I want.
//
void sample_field(
  double fields[2*D],
  double weights[8],
  double field_components[8][2*D])
{
  ASSERT_ALIGNED(fields);
  ASSERT_ALIGNED(weights);
  ASSERT_ALIGNED(field_components);
  for(int c=0; c<8; c++)
  #pragma simd
  for(int i=0;i<D*2;i++)
  {
    fields[i] += weights[c]*field_components[c][i];
  }
}

void donothing_in_serial()
{
  const Time start = get_time();
  report_time(start);
}

void donothing_in_parallel()
{
  #pragma omp parallel
  {
  const Time start = get_time();
  report_time(start);
  }
}

int main()
{
  //test_alex();
  //test_alex2();
  //test_AoS();
  printf("#\n");
  printf("# I do everything twice to separate out cache miss issues.\n");
  printf("#\n");
  printf("#\n");
  printf("# === initialization tests ===\n");
  printf("#\n");
  initialize_data();
  // do everything twice to expose cache issues
  printf("#\n");
  printf("#\n");
  printf("# === pushing tests ===\n");
  printf("#\n");
  printf("# pushing particles with no attempt at parallelization is slow:\n");
  printf("#\n");
  test_push_pcls_in_cell();
  test_push_pcls_in_cell();
  printf("#\n");
  printf("# the old SoA mover vectorizes nicely and incredibly seems to be\n");
  printf("# as fast as a simple copy even when data is in cache.\n");
  printf("# How is this possible?:\n");
  printf("#\n");
  move_bucket_old();
  move_bucket_old();
  printf("# \n");
  printf("# move particles in 8-particle transposable blocks:\n");
  move_SoA_blocks();
  move_SoA_blocks();
  printf("# \n");
  printf("# this rewrite of move_bucket_old uses three-element arrays\n");
  printf("# instead of three separate variables, which if the compiler were\n");
  printf("# intelligent would make no difference; for some bizarre reason\n");
  printf("# this takes twice as long to run when 'denom' is calculated\n");
  printf("# with double precision as for single precision, whereas in\n");
  printf("# move_bucket_old, changing to double precision increases the\n");
  printf("# cost very little.\n");
  printf("# \n");
  test_push_pcls_in_cell_SoA_vectorized();
  test_push_pcls_in_cell_SoA_vectorized();
  printf("#\n");
  printf("# the only real difference between the following and SoA\n");
  printf("# vectorization is that data must be gathered at the beginning\n");
  printf("# of the loop and scattered at the end, but unfortunately those each\n");
  printf("# seem to take as long as the push itself.\n");
  printf("#\n");
  test_push_pcls_in_cell_AoS_scatter_gather_vectorization();
  test_push_pcls_in_cell_AoS_scatter_gather_vectorization();
  printf("#\n");
  printf("# Why is this 8 times slower than SoA vectorization?  I write\n");
  printf("# loops from 0 to 7 so that the compiler would generate vector\n");
  printf("# instructions.  Is there an issue with not pipelining vector\n");
  printf("# unit data or instructions fast enough?\n");
  printf("#\n");
  test_push_pcls_in_cell_AoS_localized_vectorization();
  test_push_pcls_in_cell_AoS_localized_vectorization();
  printf("# \n");
  printf("# move particles stopping at boundary of mesh cell:\n");
  push_pcls_in_cell_SoA_stopping_at_face();
  push_pcls_in_cell_SoA_stopping_at_face();
  printf("#\n");
  printf("#\n");
  printf("# demonstrating that the time required to query the time is negligible:\n");
  printf("#\n");
  donothing_in_parallel();
  donothing_in_parallel();
  donothing_in_serial();
  donothing_in_serial();
}
