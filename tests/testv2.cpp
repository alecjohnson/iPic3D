#include <omp.h>
#include <stdio.h>
//#include "time.h" // for clock_gettime()
#include <stdint.h> // for uint64_t
#include <stdlib.h> // rand()
#include <sys/time.h>
#include <limits.h> // RAND_MAX
#include <string.h> // memcpy
#include <algorithm> // for std::max
#include <cmath> // for std::abs
// change to include this conditionally
#include <micvec.h>
#include <assert.h>
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"

#define ALIGNMENT 64
#ifdef __INTEL_COMPILER
  #define ALLOC_ALIGNED __attribute__((aligned(ALIGNMENT)));
  #define ASSUME_ALIGNED(X) __assume_aligned(X, ALIGNMENT)
#else
    #define ALLOC_ALIGNED
    #define ASSUME_ALIGNED(X)
#endif

// double-precision float
typedef double dfloat;
// float used when pushing particles
typedef dfloat pfloat;
// float used for computing reciprocal
typedef float rfloat;

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

inline dfloat get_sec(timespec start_time)
{
  timespec end_time = get_timespec();
  timespec diff_time = diff_timespec(start_time,end_time);
  return diff_time.tv_sec+1.e-9*diff_time.tv_nsec;
}

inline dfloat get_msec(timespec start_time)
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
inline dfloat time_msec()
{
  static struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec + tv.tv_usec * 1.e-6)*1.e3;
  // this is another way
  return 1.e3*omp_get_wtime();
}

// measure time using gettimeofday
//
#define Time dfloat
#define report_time(start) \
 { \
   const dfloat end = time_msec(); \
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
dfloat field_transpose[6][8] ALLOC_ALIGNED;
dfloat weights[8]ALLOC_ALIGNED;
dfloat Bx;
dfloat By;
dfloat Bz;
dfloat Ex;
dfloat Ey;
dfloat Ez;
dfloat field_components[8][8]ALLOC_ALIGNED;
dfloat fields[8]ALLOC_ALIGNED;

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
  dfloat fields[2*D],
  dfloat weights[8],
  dfloat field_components[8][2*D]);

//void sample_field(
//  dfloat fields[2*D],
//  dfloat weights[8],
//  dfloat field_components[8][2*D])__attribute__((noinline));

void test_AoS()
{
   dfloat B[4]ALLOC_ALIGNED;
   dfloat E[4]ALLOC_ALIGNED;
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

// sample from open unit interval (0,1)
dfloat sample_openu()
{
  const dfloat max_inv = 1./(dfloat(RAND_MAX)+2);
  return (dfloat(rand())+1)*max_inv;
}

// sample from closed unit interval [0,1]
dfloat usample()
{
  const dfloat RAND_MAX_inv = 1./RAND_MAX;
  return rand()*RAND_MAX_inv;
  //
}

// time step
dfloat dt;
// cell dimensions
dfloat dx;
dfloat dy;
dfloat dz;
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
  dfloat x[3];
  dfloat q;
  dfloat u[3];
  dfloat t;
  //long long ID;
 public:
  dfloat get_x(int i)const{return x[i];}
  dfloat get_u(int i)const{return u[i];}
  dfloat* fetch_x(){return x;}
  dfloat* fetch_u(){return u;}
  void set_x(dfloat* in){
    for(int i=0;i<3;i++)
      x[i]=in[i];
  }
  void set_u(dfloat* in){
    for(int i=0;i<3;i++)
      u[i]=in[i];
  }
  void set_x(int i, dfloat in){x[i]=in;}
  void set_u(int i, dfloat in){u[i]=in;}
  void init_random(int i)
  {
    for(int j=0;j<3;j++)
    {
      x[j] = usample();
      u[j] = usample();
    }
    q=usample();
    t=usample()*dt;
    //ID=i;
  }
};

// transpose each 2x2 block in the two rows of a 2x8 matrix
// used to transpose the 2x2 blocks of the 8x8 matrix
inline void trans2x2(double in1[8], double in2[8])
{
  F64vec8& data1 = (F64vec8&)in1[0];
  F64vec8& data2 = (F64vec8&)in2[0];
  // copy the data that we will first overwrite
  const F64vec8 buff1 = data1;

  // data1: copy odds from evens of data2,
  //   i.e.,
  // for(int i=1; i<8; i+=2) data1[i] = data2[i-1];
  //
  const __mmask16 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
  // replace unmasked data1 with swizzled data2
  // (hopefully compiler will see that this can be done in
  // a single vector instruction)
  data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
    rmask_01010101, __m512d(data2),_MM_SWIZ_REG_CDAB));
  //cout << "swizzle for pattern 'cdab' with rmask=01010101\n"
  //  << data1 << endl;

  // data2: copy evens from odds of buff1,
  //   i.e.,
  // for(int i=0; i<7; i+=2) data2[i] = buff1[i+1];
  //
  const __mmask16 rmask_10101010 = _mm512_int2mask(0x55); // mask=01010101
  // replace unmasked data2 with swizzled buff1
  data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
    rmask_10101010, __m512d(buff1),_MM_SWIZ_REG_CDAB));
  //cout << "swizzle for pattern 'cdab' with rmask=10101010\n"
  //  << data2 << endl;
}

// transpose 1x2 elements of each 2x4 block in the two rows of a 2x8 matrix
// used to transpose the 2x2 elements of the 4x4 blocks of the 8x8 matrix
inline void trans4x4(double in1[8], double in2[8])
{
  F64vec8& data1 = (F64vec8&)in1[0];
  F64vec8& data2 = (F64vec8&)in2[0];
  // copy the data that we will first overwrite
  const F64vec8 buff1 = data1;

  // data1: copy odd 1x2 elements from even 1x2 elements of data2,
  //   i.e.,
  // for(int i=2; i<8; i++) if((i/2)%2) data1[i] = data2[i-2];
  //
  // replace unmasked data1 with swizzled data2
  const __mmask16 rmask_00110011 = _mm512_int2mask(0xcc); // mask = 11001100
  data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
    rmask_00110011, __m512d(data2),_MM_SWIZ_REG_BADC));
  //cout << "swizzle for pattern 'badc' with rmask=00110011\n"
  //  << outp1 << endl;

  // data2: copy even 1x2 elements from odd 1x2 from odds of buff1,
  //   i.e.,
  // for(int i=0; i<6; i++) if(!((i/2)%2)) data2[i] = buff1[i+2];
  //
  const __mmask16 rmask_11001100 = _mm512_int2mask(0x33); // mask=00110011
  data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
    rmask_11001100, __m512d(buff1),_MM_SWIZ_REG_BADC));
  //cout << "swizzle for pattern 'badc' with rmask=11001100\n"
  //  << outp2 << endl;
}

// transpose 1x4 elements of each 2x4 block in the two rows of a 2x8 matrix
// used to transpose the 4x4 elements of the 8x8 matrix
void trans8x8(double in1[8], double in2[8])
{
  __m512& data1 = (__m512&)in1[0];
  __m512& data2 = (__m512&)in2[0];
  // copy the data that we will first overwrite
  const __m512 buff1 = data1;

  // for swizzle intel supports 256-bit lanes with 64-bit elements
  // (as well as 128-bit lanes with 32-bit elements),
  // but for shuffle intel only supports 128-bit lanes with 32-bit
  // elements, so we must cast to 32-bit data types, do the
  // shuffle, and cast back.

  // data1: copy low (odd) 1x4 element from even (high) 1x4 element of data2,
  //   i.e.,
  // for(int i=0; i<4; i++) data1[i+4] = data2[i];
  //
  // replace unmasked data1 with shuffled data2
  const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
  data1 = _mm512_mask_permute4f128_ps(data1, rmask_00001111, data2,_MM_PERM_BADC);

  // data2: copy high (even) 1x4 element from odd (low) 1x4 element of buff1,
  //   i.e.,
  // for(int i=0; i<4; i++) data2[i] = buff1[i+4];
  //
  // replace unmasked data2 with shuffled data1
  const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
  data2 = _mm512_mask_permute4f128_ps(data2, rmask_11110000, buff1,_MM_PERM_BADC);
}

// transpose with blocked algorithm in
// 24=8*3 = 8*log_2(8) vector instructions
// (ignoring 12=4*3 copies made to
// buffers to reduce number of registers needed)
//
inline void transpose_8x8_double(double data[8][8])
{
  ASSUME_ALIGNED(data);
  // 1. transpose each 2x2 block.
  for(int i=0; i<8; i+=2)
    trans2x2(data[i], data[i+1]);
  // 2. transpose each 4x4 block of 2x2 elements
  trans4x4(data[0], data[2]);
  trans4x4(data[1], data[3]);
  trans4x4(data[4], data[6]);
  trans4x4(data[5], data[7]);
  // 3. swap lower left and upper right 4x4 elements
  for(int i=0; i<4; i+=1)
    trans8x8(data[i], data[i+4]);
}

// 8 particles
class PclBlock
{
  dfloat data[8][8];
 public:
  const void* get_data_memory()const
  {
    return &data[0][0];
  }
  void naive_transpose()
  {
    dfloat temp[8][8]ALLOC_ALIGNED;
    for(int i=0;i<8;i++)
    for(int j=0;j<8;j++)
    {
      temp[i][j] = data[j][i];
    }
    memcpy(&data[0][0],&temp[0][0],sizeof(dfloat)*8*8);
  }
  // convert between AoS and SoA
  void transpose()
  {
    //naive_transpose();
    transpose_8x8_double(data);
  }
  // assumes AoS
  SpeciesPcl& fetch_pcl(int p){return (SpeciesPcl&)data[p];}
  // assumes SoA
  dfloat* fetch_x(){return data[0];}
  dfloat* fetch_y(){return data[1];}
  dfloat* fetch_z(){return data[2];}
  dfloat* fetch_q(){return data[3];}
  dfloat* fetch_u(){return data[4];}
  dfloat* fetch_v(){return data[5];}
  dfloat* fetch_w(){return data[6];}
  dfloat* fetch_t(){return data[7];}
  // we will not track IDs
  //long long* fetch_ID(){return (long long*) data[7];}
};
bool operator== (const PclBlock &lhs, const PclBlock &rhs)
{
  const void* lhsdata = lhs.get_data_memory();
  const void* rhsdata = rhs.get_data_memory();
  return !memcmp(lhsdata, rhsdata, sizeof(dfloat)*8*8);
  //for(int i=0; i<8; i++)
  //for(int j=0; j<8; j++)
  //  lhs.data[i][j]==rhs.data[i][j];
}

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

const int NUMPCLS=NPB*256;
const int NUMBLKS=(NUMPCLS-1)/8+1;
SpeciesPcl pcls[NUMPCLS]ALLOC_ALIGNED;
SpeciesPcl pcls2[NUMPCLS]ALLOC_ALIGNED;
PclBlock pclBlocks[NUMBLKS]ALLOC_ALIGNED;
dfloat x[NUMPCLS]ALLOC_ALIGNED;
dfloat y[NUMPCLS]ALLOC_ALIGNED;
dfloat z[NUMPCLS]ALLOC_ALIGNED;
dfloat u[NUMPCLS]ALLOC_ALIGNED;
dfloat v[NUMPCLS]ALLOC_ALIGNED;
dfloat w[NUMPCLS]ALLOC_ALIGNED;
dfloat q[NUMPCLS]ALLOC_ALIGNED;
long long ID[NUMPCLS]ALLOC_ALIGNED;
// portion of dt remaining for particle
dfloat t[NUMPCLS]ALLOC_ALIGNED;
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
     pcls[i].t = t[i];
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
    t[i] = pcls[i].t;
  }
  report_time(start);
  }
}

// Sebastian's revised version
//
void copy_AoS_particles_to_another_array_rev()
{
  uint64_t *p1 = (uint64_t *) pcls;  ASSUME_ALIGNED(p1);
  uint64_t *p2 = (uint64_t *) pcls2; ASSUME_ALIGNED(p2);
  #pragma omp parallel
  {
  //printf("sizeof(SpeciesPcl): %d\n", sizeof(SpeciesPcl));
  const Time start = get_time();
  // With schedule() we make sure that every thread gets aligned
  // consecutive chunk of pcls to copy
  #pragma omp for schedule(static,NUMPCLS*4) nowait
  for(int p=0;p<NUMPCLS*8;p++)
  {
    //ASSUME_ALIGNED(&pcls2[0]);
    //ASSUME_ALIGNED(&pcls[0]);
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
    ASSUME_ALIGNED(&pcls2[0]);
    ASSUME_ALIGNED(&pcls[0]);
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
      ASSUME_ALIGNED(pcl);
      ASSUME_ALIGNED(&pcls[0]);
      (*pcl) = pcls[b*8+p];
    }
    PclBlock tmpBlock1 = pclBlock;
    PclBlock tmpBlock2 = pclBlock;
    tmpBlock1.naive_transpose();
    tmpBlock2.transpose();
    assert(tmpBlock1==tmpBlock2);
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
  dt = sample_openu();
  dprint(dt);
  dx = usample();
  dy = usample();
  dz = usample();
  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].init_random(i);
    // ensure that t is no greater than dt
    t[i] = dt*usample();
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
  dfloat dto2 = usample();
  dfloat qdto2mc = usample();
  dfloat cellstart[D];
  dfloat dx_inv[D];
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
    ASSUME_ALIGNED(pcl);
    dfloat xorig[D]ALLOC_ALIGNED;
    dfloat xavg[D]ALLOC_ALIGNED;
    dfloat uorig[D]ALLOC_ALIGNED;

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
    dfloat B[D]ALLOC_ALIGNED;
    dfloat E[D]ALLOC_ALIGNED;
    {
      dfloat fields[D*2]ALLOC_ALIGNED;
      dfloat weights[8]ALLOC_ALIGNED;
      {
        dfloat w[2][D]ALLOC_ALIGNED;
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
    dfloat uavg[D]ALLOC_ALIGNED;
    {
      dfloat Om[D]ALLOC_ALIGNED;
      dfloat denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      {
        dfloat omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1/rfloat(omsq_p1);
      }
      dfloat ut[D]ALLOC_ALIGNED;
      dfloat udotOm;
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
  dfloat dto2 = usample();
  dfloat qdto2mc = usample();
  dfloat cellstart[D];
  dfloat dx_inv[D];
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
    dfloat xorig[D];
    dfloat xavg[D];
    dfloat uorig[D];

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
    dfloat B[D];
    dfloat E[D];
    {
      dfloat fields[D*2];
      dfloat weights[8];
      {
        dfloat w[2][D];
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
    dfloat uavg[D];
    {
      dfloat Om[D];
      dfloat denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      //const dfloat omsq_p1 = 1.0+(Om[0] * Om[0] + Om[1] * Om[1] + Om[2] * Om[2]);
      //const dfloat denom = 1/rfloat(omsq_p1);
      {
        dfloat omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1/rfloat(omsq_p1);
        //denom = 1.0/omsq_p1;
      }
      dfloat ut[D];
      dfloat udotOm;
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
// calculated with dfloat precision.
//
void move_bucket_old()
{
    const dfloat dto2 = usample();
    const dfloat qdto2mc = usample();
    const dfloat cellstartx=usample();
    const dfloat cellstarty=usample();
    const dfloat cellstartz=usample();
    const dfloat dx_inv=usample();
    const dfloat dy_inv=usample();
    const dfloat dz_inv=usample();
    #pragma omp parallel 
    {
        const Time start = get_time();
        #pragma omp for
        #pragma simd
        for(int pidx = 0; pidx < NUMPCLS; pidx++)
        {
            // copy the particle
            const dfloat xorig = x[pidx];
            const dfloat yorig = y[pidx];
            const dfloat zorig = z[pidx];
            const dfloat uorig = u[pidx];
            const dfloat vorig = v[pidx];
            const dfloat worig = w[pidx];

            // initialize xavg to xorig
            dfloat xavg = x[pidx];
            dfloat yavg = y[pidx];
            dfloat zavg = z[pidx];

            // compute weights for field components
            //
            dfloat weights[8];
            // fraction of the distance from the left of the cell
            const dfloat w0x = dx_inv*(xavg - cellstartx);
            const dfloat w0y = dy_inv*(yavg - cellstarty);
            const dfloat w0z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the right
            const dfloat w1x = 1.-w0x;
            const dfloat w1y = 1.-w0y;
            const dfloat w1z = 1.-w0z;
            const dfloat weight00 = w0x*w0y;
            const dfloat weight01 = w0x*w1y;
            const dfloat weight10 = w1x*w0y;
            const dfloat weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            // Interpolate the electromagnetic field
            //
            dfloat Exl=0.0, Eyl=0.0, Ezl=0.0;
            dfloat Bxl=0.0, Byl=0.0, Bzl=0.0;
            for(int c=0; c<8; c++)
            {
                Bxl += weights[c] * field_components[c][0];
                Byl += weights[c] * field_components[c][1];
                Bzl += weights[c] * field_components[c][2];
                Exl += weights[c] * field_components[c][3];
                Eyl += weights[c] * field_components[c][4];
                Ezl += weights[c] * field_components[c][5];
            }

            const dfloat Omx = qdto2mc*Bxl;
            const dfloat Omy = qdto2mc*Byl;
            const dfloat Omz = qdto2mc*Bzl;

            // end interpolation
            const dfloat omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const dfloat denom = 1/rfloat(omsq_p1);
            // solve the position equation
            const dfloat ut = uorig + qdto2mc * Exl;
            const dfloat vt = vorig + qdto2mc * Eyl;
            const dfloat wt = worig + qdto2mc * Ezl;
            //const dfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const dfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const dfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const dfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const dfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
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
static inline pfloat get_truncation_ratio(
  //int&direction,
  pfloat dXwant, pfloat dXfull,
  pfloat dYwant, pfloat dYfull,
  pfloat dZwant, pfloat dZfull)
{
  const pfloat motion_freedom = 1.e-2;
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

  const pfloat denominator = dXfull*dYfull*dZfull;
  const pfloat dXprod = dXwant*dYfull*dZfull;
  const pfloat dYprod = dYwant*dXfull*dZfull;
  const pfloat dZprod = dZwant*dXfull*dYfull;
  pfloat numerator; // = denominator;
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
  return numerator/rfloat(denominator);
}

// in which direction are we most outside?
// +/-1: +/-X
// +/-2: +/-Y
// +/-4: +/-Z
static inline int get_direction(pfloat Xpos, pfloat Ypos, pfloat Zpos)
{
  int direction;
  const pfloat aX = std::abs(Xpos);
  const pfloat aY = std::abs(Ypos);
  const pfloat aZ = std::abs(Zpos);
  //
  // This way requires 8 comparisons and 7 masked assignments
  //
  //if(aX > aY) direction = (aX>aZ) ? 1 : 4;
  //else direction = (aY>aZ) ? 2 : 4;
  //const pfloat low = max(xpos,ypos,zpos);
  //const pfloat hgh = min(xpos,ypos,zpos);
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

void push_pcls_in_cell_SoA_stopping_at_face()
{
    // time step resolution (analogous to FLT_MIN,
    // defines a limit on precision)
    const pfloat dt_min = 1e-6*dt;
    // shortest allowed subcycle time step
    // (assumed to be no less than dt_min)
    const pfloat min_dt = 1e-5*dt;
    //const dfloat dto2 = usample();
    //const dfloat qdto2mc = usample();
    const dfloat qo2mc = usample();
    const dfloat dx_over_two = dx/2;
    const dfloat dy_over_two = dy/2;
    const dfloat dz_over_two = dz/2;
    const dfloat two_over_dx = 2/dx;
    const dfloat two_over_dy = 2/dy;
    const dfloat two_over_dz = 2/dz;
    const dfloat xmiddle = usample(); // position of middle of cell
    const dfloat ymiddle = usample(); // position of middle of cell
    const dfloat zmiddle = usample(); // position of middle of cell
    #pragma omp parallel 
    {
        const Time start = get_time();
        #pragma omp for
        #pragma simd // vectorlength(16) // why doesn't this help?
        for(int pidx = 0; pidx < NUMPCLS; pidx++)
        {
          // copy the particle
          //
          // x is physical position
          // X is position is in canonical coordinates (-1 <=~ X <=~ 1)
          //
          //const pfloat Xorig = X[pidx];
          //const pfloat Yorig = Y[pidx];
          //const pfloat Zorig = Z[pidx];
          const pfloat Xorig = (x[pidx]-xmiddle)*two_over_dx;
          const pfloat Yorig = (y[pidx]-ymiddle)*two_over_dy;
          const pfloat Zorig = (z[pidx]-zmiddle)*two_over_dz;
          // u is physical velocity
          const pfloat uorig = u[pidx];
          const pfloat vorig = v[pidx];
          const pfloat worig = w[pidx];
          const pfloat torig = t[pidx];
          //assert_le(torig, dt);

          // The computed time step needs to be dfloat
          // precision up to the point where the calculation is
          // unique for every particle.
          //
          // compute time remaining for particle until
          // next synchonization point.  Note that if torig==0
          // then there is no loss of precision at this point.
          dfloat dtpcl = dt-torig;
          // initialize subcycle time to be remaining time
          // (used in first iteration of iterative solver)
          //
          dfloat dtcycle = dtpcl;

          // initialize xavg to xorig
          pfloat Xavg = Xorig;
          pfloat Yavg = Yorig;
          pfloat Zavg = Zorig;

          // purpose of iterative solver is to find
          // dtcycle, Xavg.., and uavg...
          pfloat uavg, vavg, wavg;

          // this is the part that must vectorize
          // #pragma omp simd
          for(int niter=0;niter<NiterMover;niter++)
          {
            // compute weights for field components
            //
            pfloat weights[8];
            // fraction of the distance from the left of the cell
            const pfloat w0x = 0.5*Xavg + 0.5;
            const pfloat w0y = 0.5*Yavg + 0.5;
            const pfloat w0z = 0.5*Zavg + 0.5;
            // fraction of distance from the right
            const dfloat w1x = 1.-w0x;
            const dfloat w1y = 1.-w0y;
            const dfloat w1z = 1.-w0z;
            const dfloat weight00 = w0x*w0y;
            const dfloat weight01 = w0x*w1y;
            const dfloat weight10 = w1x*w0y;
            const dfloat weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            pfloat Exl = 0.0;
            pfloat Eyl = 0.0;
            pfloat Ezl = 0.0;
            pfloat Bxl = 0.0;
            pfloat Byl = 0.0;
            pfloat Bzl = 0.0;

            // field_components is dfloat precision; loss of
            // precision at this point is expected to mitigated
            // by the large number of particles.  When we sum
            // moments we will go back to dfloat precision.
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

            const dfloat qdto2mc = qo2mc*dtcycle;
            const pfloat Omx = qdto2mc*Bxl;
            const pfloat Omy = qdto2mc*Byl;
            const pfloat Omz = qdto2mc*Bzl;

            // end interpolation
            const pfloat omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const pfloat denom = 1/rfloat(omsq_p1);
            // solve the position equation
            const pfloat ut = uorig + qdto2mc * Exl;
            const pfloat vt = vorig + qdto2mc * Eyl;
            const pfloat wt = worig + qdto2mc * Ezl;
            //const dfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const pfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const pfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const pfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

            // stop the particle at the cell boundary
            //
            // compute the displacement assuming the particle is not stopped
            //
            const pfloat dxpcl = dtcycle*uavg;
            const pfloat dypcl = dtcycle*vavg;
            const pfloat dzpcl = dtcycle*wavg;
            const pfloat dXpcl = dxpcl*two_over_dx;
            const pfloat dYpcl = dypcl*two_over_dy;
            const pfloat dZpcl = dzpcl*two_over_dz;
            const pfloat Xnew = Xorig + dXpcl;
            const pfloat Ynew = Yorig + dYpcl;
            const pfloat Znew = Zorig + dZpcl;
            //
            // compute the factor by which the motion
            // (time step) must be multiplied for the particle
            // to stop at the cell boundary.
            //
            // 1. compute the distance moved.
            //
            const pfloat dXmag = std::abs(dXpcl);
            const pfloat dYmag = std::abs(dYpcl);
            const pfloat dZmag = std::abs(dZpcl);
            //
            // 2. compute the distance to the wall
            //    (if moving away then allow no motion)
            //
            pfloat dXwall, dYwall, dZwall;
            const pfloat hghXwall=1., lowXwall=-1.;
            const pfloat hghYwall=1., lowYwall=-1.;
            const pfloat hghZwall=1., lowZwall=-1.;
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
            // dfloat the distance and therefore should not cause
            // a problem in the final result; a subsequent iteration
            // should then be able to take the particle to the wall...).
            //dXwall = dXmag + 1 - abs(Xnew);
            //dXwall = max(abs(Xorig)-1, dXwall);
            //
            // 3. compute the ratio to truncate motion
            //
            const pfloat ratio = get_truncation_ratio(
              dXwall, dXmag,
              dYwall, dYmag,
              dZwall, dZmag);
            //
            // 4. truncate or adjust the time step and motion accordingly
            //
            const pfloat dtproposed = dtcycle*ratio;
            // enforce a minimum dtcycle
            dtcycle = std::max(min_dt, dtproposed);
            // but cap final time at synchronization time
            // (at which point dfloat precision is here restored)
            dtcycle = std::min(dtpcl, dtcycle);

            // apply the corrected time step
            Xavg = Xorig + 0.5*dtcycle*uavg;
            Yavg = Yorig + 0.5*dtcycle*vavg;
            Zavg = Zorig + 0.5*dtcycle*wavg;
          }

          // update particle after the last iteration
          const pfloat Xend = 2.0 * Xavg - Xorig;
          const pfloat Yend = 2.0 * Yavg - Yorig;
          const pfloat Zend = 2.0 * Zavg - Zorig;
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
            t[pidx] += dtcycle;
          }

          // if the particle is done being moved then it
          // can be presumed inside the box
          if(t[pidx] > (dt-dt_min))
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

// try to take a time step that stops particle at a face
//
void push_SoA_blocks_stopping_at_face()
{
    // time step resolution (analogous to FLT_MIN,
    // defines a limit on precision)
    const pfloat dt_min = 1e-6*dt;
    // shortest allowed subcycle time step
    // (assumed to be no less than dt_min)
    const pfloat min_dt = 1e-5*dt;
    //const dfloat dto2 = usample();
    //const dfloat qdto2mc = usample();
    const dfloat qo2mc = usample();
    const dfloat dx_over_two = dx/2;
    const dfloat dy_over_two = dy/2;
    const dfloat dz_over_two = dz/2;
    const dfloat two_over_dx = 2/dx;
    const dfloat two_over_dy = 2/dy;
    const dfloat two_over_dz = 2/dz;
    const dfloat xmiddle = usample(); // position of middle of cell
    const dfloat ymiddle = usample(); // position of middle of cell
    const dfloat zmiddle = usample(); // position of middle of cell
    #pragma omp parallel 
    {
      const Time start = get_time();
      #pragma omp for
      for(int bidx = 0; bidx < NUMBLKS; bidx++)
      {
        PclBlock& pclBlock = pclBlocks[bidx];
        pclBlock.transpose();
        dfloat* x = pclBlock.fetch_x();
        dfloat* y = pclBlock.fetch_y();
        dfloat* z = pclBlock.fetch_z();
        dfloat* u = pclBlock.fetch_u();
        dfloat* v = pclBlock.fetch_v();
        dfloat* w = pclBlock.fetch_w();
        dfloat* t = pclBlock.fetch_t();
        ASSUME_ALIGNED(x);
        ASSUME_ALIGNED(y);
        ASSUME_ALIGNED(z);
        ASSUME_ALIGNED(u);
        ASSUME_ALIGNED(v);
        ASSUME_ALIGNED(w);
        ASSUME_ALIGNED(t);

        // push all the particles in the block
        #pragma simd
        for(int pidx = 0; pidx < 8; pidx++)
        {
          // copy the particle
          //
          // x is physical position
          // X is position is in canonical coordinates (-1 <=~ X <=~ 1)
          //
          //const pfloat Xorig = X[pidx];
          //const pfloat Yorig = Y[pidx];
          //const pfloat Zorig = Z[pidx];
          const pfloat Xorig = (x[pidx]-xmiddle)*two_over_dx;
          const pfloat Yorig = (y[pidx]-ymiddle)*two_over_dy;
          const pfloat Zorig = (z[pidx]-zmiddle)*two_over_dz;
          // u is physical velocity
          const pfloat uorig = u[pidx];
          const pfloat vorig = v[pidx];
          const pfloat worig = w[pidx];
          const pfloat torig = t[pidx];
          //assert_le(torig, dt);

          // The computed time step needs to be dfloat
          // precision up to the point where the calculation is
          // unique for every particle.
          //
          // compute time remaining for particle until
          // next synchonization point.  Note that if torig==0
          // then there is no loss of precision at this point.
          dfloat dtpcl = dt-torig;
          // initialize subcycle time to be remaining time
          // (used in first iteration of iterative solver)
          //
          dfloat dtcycle = dtpcl;

          // initialize xavg to xorig
          pfloat Xavg = Xorig;
          pfloat Yavg = Yorig;
          pfloat Zavg = Zorig;

          // purpose of iterative solver is to find
          // dtcycle, Xavg.., and uavg...
          pfloat uavg, vavg, wavg;

          // this is the part that must vectorize
          // #pragma omp simd
          for(int niter=0;niter<NiterMover;niter++)
          {
            // compute weights for field components
            //
            pfloat weights[8];
            // fraction of the distance from the left of the cell
            const pfloat w0x = 0.5*Xavg + 0.5;
            const pfloat w0y = 0.5*Yavg + 0.5;
            const pfloat w0z = 0.5*Zavg + 0.5;
            // fraction of distance from the right
            const dfloat w1x = 1.-w0x;
            const dfloat w1y = 1.-w0y;
            const dfloat w1z = 1.-w0z;
            const dfloat weight00 = w0x*w0y;
            const dfloat weight01 = w0x*w1y;
            const dfloat weight10 = w1x*w0y;
            const dfloat weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            pfloat Exl = 0.0;
            pfloat Eyl = 0.0;
            pfloat Ezl = 0.0;
            pfloat Bxl = 0.0;
            pfloat Byl = 0.0;
            pfloat Bzl = 0.0;

            // field_components is dfloat precision; loss of
            // precision at this point is expected to mitigated
            // by the large number of particles.  When we sum
            // moments we will go back to dfloat precision.
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

            const dfloat qdto2mc = qo2mc*dtcycle;
            const pfloat Omx = qdto2mc*Bxl;
            const pfloat Omy = qdto2mc*Byl;
            const pfloat Omz = qdto2mc*Bzl;

            // end interpolation
            const pfloat omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const pfloat denom = 1/rfloat(omsq_p1);
            // solve the position equation
            const pfloat ut = uorig + qdto2mc * Exl;
            const pfloat vt = vorig + qdto2mc * Eyl;
            const pfloat wt = worig + qdto2mc * Ezl;
            //const dfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const pfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const pfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const pfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

            // stop the particle at the cell boundary
            //
            // compute the displacement assuming the particle is not stopped
            //
            const pfloat dxpcl = dtcycle*uavg;
            const pfloat dypcl = dtcycle*vavg;
            const pfloat dzpcl = dtcycle*wavg;
            const pfloat dXpcl = dxpcl*two_over_dx;
            const pfloat dYpcl = dypcl*two_over_dy;
            const pfloat dZpcl = dzpcl*two_over_dz;
            const pfloat Xnew = Xorig + dXpcl;
            const pfloat Ynew = Yorig + dYpcl;
            const pfloat Znew = Zorig + dZpcl;
            //
            // compute the factor by which the motion
            // (time step) must be multiplied for the particle
            // to stop at the cell boundary.
            //
            // 1. compute the distance moved.
            //
            const pfloat dXmag = std::abs(dXpcl);
            const pfloat dYmag = std::abs(dYpcl);
            const pfloat dZmag = std::abs(dZpcl);
            //
            // 2. compute the distance to the wall
            //    (if moving away then allow no motion)
            //
            pfloat dXwall, dYwall, dZwall;
            const pfloat hghXwall=1., lowXwall=-1.;
            const pfloat hghYwall=1., lowYwall=-1.;
            const pfloat hghZwall=1., lowZwall=-1.;
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
            // dfloat the distance and therefore should not cause
            // a problem in the final result; a subsequent iteration
            // should then be able to take the particle to the wall...).
            //dXwall = dXmag + 1 - abs(Xnew);
            //dXwall = max(abs(Xorig)-1, dXwall);
            //
            // 3. compute the ratio to truncate motion
            //
            const pfloat ratio = get_truncation_ratio(
              dXwall, dXmag,
              dYwall, dYmag,
              dZwall, dZmag);
            //
            // 4. truncate or adjust the time step and motion accordingly
            //
            const pfloat dtproposed = dtcycle*ratio;
            // enforce a minimum dtcycle
            dtcycle = std::max(min_dt, dtproposed);
            // but cap final time at synchronization time
            // (at which point dfloat precision is here restored)
            dtcycle = std::min(dtpcl, dtcycle);

            // apply the corrected time step
            Xavg = Xorig + 0.5*dtcycle*uavg;
            Yavg = Yorig + 0.5*dtcycle*vavg;
            Zavg = Zorig + 0.5*dtcycle*wavg;
          }

          // update particle after the last iteration
          const pfloat Xend = 2.0 * Xavg - Xorig;
          const pfloat Yend = 2.0 * Yavg - Yorig;
          const pfloat Zend = 2.0 * Zavg - Zorig;
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
            t[pidx] += dtcycle;
          }

          // if the particle is done being moved then it
          // can be presumed inside the box
          if(t[pidx] > (dt-dt_min))
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
        pclBlock.transpose();
      }
      report_time(start);
    }
}

// move particles in 8-particle transposable blocks
void move_SoA_blocks()
{
    const dfloat dto2 = usample();
    const dfloat qdto2mc = usample();
    const dfloat cellstartx=usample();
    const dfloat cellstarty=usample();
    const dfloat cellstartz=usample();
    const dfloat dx_inv=usample();
    const dfloat dy_inv=usample();
    const dfloat dz_inv=usample();
    #pragma omp parallel 
    {
      const Time start = get_time();
      #pragma omp for
      for(int bidx = 0; bidx < NUMBLKS; bidx++)
      {
        PclBlock& pclBlock = pclBlocks[bidx];
        pclBlock.transpose();
        dfloat* x = pclBlock.fetch_x();
        dfloat* y = pclBlock.fetch_y();
        dfloat* z = pclBlock.fetch_z();
        dfloat* u = pclBlock.fetch_u();
        dfloat* v = pclBlock.fetch_v();
        dfloat* w = pclBlock.fetch_w();
        ASSUME_ALIGNED(x);
        ASSUME_ALIGNED(y);
        ASSUME_ALIGNED(z);
        ASSUME_ALIGNED(u);
        ASSUME_ALIGNED(v);
        ASSUME_ALIGNED(w);
        #pragma omp simd
        for(int pidx = 0; pidx < 8; pidx++)
        {
            // copy the particle
            const dfloat xorig = x[pidx];
            const dfloat yorig = y[pidx];
            const dfloat zorig = z[pidx];
            const dfloat uorig = u[pidx];
            const dfloat vorig = v[pidx];
            const dfloat worig = w[pidx];

            // initialize xavg to xorig
            dfloat xavg = x[pidx];
            dfloat yavg = y[pidx];
            dfloat zavg = z[pidx];

            // compute weights for field components
            //
            dfloat weights[8];
            // fraction of the distance from the left of the cell
            const dfloat w0x = dx_inv*(xavg - cellstartx);
            const dfloat w0y = dy_inv*(yavg - cellstarty);
            const dfloat w0z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the right
            const dfloat w1x = 1-w0x;
            const dfloat w1y = 1-w0y;
            const dfloat w1z = 1-w0z;
            const dfloat weight00 = w0x*w0y;
            const dfloat weight01 = w0x*w1y;
            const dfloat weight10 = w1x*w0y;
            const dfloat weight11 = w1x*w1y;
            weights[0] = weight00*w0z; // weight000
            weights[1] = weight00*w1z; // weight001
            weights[2] = weight01*w0z; // weight010
            weights[3] = weight01*w1z; // weight011
            weights[4] = weight10*w0z; // weight100
            weights[5] = weight10*w1z; // weight101
            weights[6] = weight11*w0z; // weight110
            weights[7] = weight11*w1z; // weight111

            dfloat Exl = 0.0;
            dfloat Eyl = 0.0;
            dfloat Ezl = 0.0;
            dfloat Bxl = 0.0;
            dfloat Byl = 0.0;
            dfloat Bzl = 0.0;

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

            const dfloat Omx = qdto2mc*Bxl;
            const dfloat Omy = qdto2mc*Byl;
            const dfloat Omz = qdto2mc*Bzl;

            // end interpolation
            const dfloat omsq_p1 = 1.0 + (Omx * Omx + Omy * Omy + Omz * Omz);
            const dfloat denom = 1/rfloat(omsq_p1);
            // solve the position equation
            const dfloat ut = uorig + qdto2mc * Exl;
            const dfloat vt = vorig + qdto2mc * Eyl;
            const dfloat wt = worig + qdto2mc * Ezl;
            //const dfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
            const dfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
            // solve the velocity equation
            const dfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
            const dfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
            const dfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
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
  dfloat dto2 = usample();
  dfloat qdto2mc = usample();
  dfloat cellstart[D];
  dfloat dx_inv[D];
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
    ASSUME_ALIGNED(pcl);
    dfloat xorig[D]ALLOC_ALIGNED;
    dfloat xavg[D]ALLOC_ALIGNED;
    dfloat uorig[D]ALLOC_ALIGNED;

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
    dfloat B[D]ALLOC_ALIGNED;
    dfloat E[D]ALLOC_ALIGNED;
    {
      dfloat fields[D*2]ALLOC_ALIGNED;
      dfloat weights[8]ALLOC_ALIGNED;
      if(true)
      {
        const dfloat w1x = dx_inv[0]*(xavg[0] - cellstart[0]);
        const dfloat w1y = dx_inv[1]*(xavg[1] - cellstart[1]);
        const dfloat w1z = dx_inv[2]*(xavg[2] - cellstart[2]);
        const dfloat w0x = 1-w1x;
        const dfloat w0y = 1-w1y;
        const dfloat w0z = 1-w1z;
        const dfloat weight00 = w0x*w0y;
        const dfloat weight01 = w0x*w1y;
        const dfloat weight10 = w1x*w0y;
        const dfloat weight11 = w1x*w1y;
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
        dfloat w[2][D]ALLOC_ALIGNED;
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
    dfloat uavg[D]ALLOC_ALIGNED;
    {
      dfloat Om[D]ALLOC_ALIGNED;
      dfloat denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      {
        dfloat omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1/rfloat(omsq_p1);
      }
      dfloat ut[D]ALLOC_ALIGNED;
      dfloat udotOm;
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

// return concatenation of low half of lo and low half of hi
inline __m512 cat_low_halves(__m512 lo, __m512 hi)
{
  const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
  // low half uses lo, and high half uses low half of hi
  return _mm512_mask_permute4f128_ps(lo, rmask_00001111, hi,_MM_PERM_BADC);
}

// return concatenation of hgh half of lo and hgh half of hi
inline __m512 cat_hgh_halves(__m512 lo, __m512 hi)
{
  const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
  // high half uses hi, and low half uses high half of lo
  return _mm512_mask_permute4f128_ps(hi, rmask_11110000, lo,_MM_PERM_BADC);
}

inline F64vec8 cat_low_halves(void* lo, void* hi)
{ return (F64vec8) cat_low_halves(*(__m512*)lo, *(__m512*)hi); }
inline F64vec8 cat_hgh_halves(void* lo, void* hi)
{ return (F64vec8) cat_hgh_halves(*(__m512*)lo, *(__m512*)hi); }
//
inline F64vec8 cat_low_halves(F64vec8 lo, F64vec8 hi)
{ return (F64vec8) cat_low_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }
inline F64vec8 cat_hgh_halves(F64vec8 lo, F64vec8 hi)
{ return (F64vec8) cat_hgh_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }

// copy low vector of src to low vector of dst
inline void copy012to012(F64vec8& dst, F64vec8 src)
{
  const __mmask8 rmask_11100000 = _mm512_int2mask(0x07); // mask=00000111
  // masked copy
  dst = _mm512_mask_mov_pd(dst,rmask_11100000,src);
}
// copy hgh vector of src to hgh vector of dst
inline void copy456to456(F64vec8& dst, F64vec8 src)
{
  const __mmask8 rmask_00001110 = _mm512_int2mask(0x70); // mask=01110000
  dst = _mm512_mask_mov_pd(dst,rmask_00001110,src);
}
// copy low vector of src to hgh vector of dst
inline void copy012to456(F64vec8& dst, F64vec8 src)
{
  // mask=0000000000111111 because 32-bit elements
  const __mmask16 rmask_1111110000000000 = _mm512_int2mask(0x03f);
  // write the first three elements of src beginning at dst[4]
  _mm512_mask_packstorelo_epi32(&dst[4],rmask_1111110000000000,*(__m512i*)&src);
}
// copy hgh vector of src to low vector of dst
inline void copy456to012(F64vec8& dst, F64vec8 src)
{
  // mask=0011111111111111 because 32-bit elements
  const __mmask16 rmask_1111111111111100 = _mm512_int2mask(0x3fff);
  // write all but the highest element to the left of dst[4].
  _mm512_mask_packstorehi_epi32(&dst[4],rmask_1111111111111100,*(__m512i*)&src);
}

// copy vectors of src to vectors of dst
inline void copy012and456(F64vec8& dst, F64vec8 src)
{
  const __mmask8 rmask_11101110 = _mm512_int2mask(0x77); //mask=01110111
  dst = _mm512_mask_mov_pd(dst,rmask_11101110,src);
}

// broadcast s0 into mask=0 slots, s1 into mask=1 slots
inline F64vec8 broadcast_mask_blend(double s0, double s1, __mmask8 mask)
{
  return _mm512_mask_blend_pd(mask, F64vec8(s0), F64vec8(s1));
  // broadcast s0 into vector
  //const F64vec8 t1 = _mm512_extload_pd(&s0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
  // broadcast s1 into same vector with mask
  //return (F64vec8) _mm512_mask_extload_pd(t1, mask, &s1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
}

inline F64vec8 cross_product(F64vec8 u, F64vec8 v)
{
  F64vec8 us = u.dacb();
  F64vec8 vs = v.dacb();
  return us*vs.dacb() - us.dacb()*vs;
}

inline F64vec8 reciprocal_F64vec8(F64vec8 arg)
{
  // is there a faster way?
  return F64vec8(1.)/arg;
  // convert to single precision for faster reciprocal
  //
  //F32vec16 arg32 = _mm512_cvtpd_pslo(arg); 
  //F32vec16 arg32inv = F32vec16(1.)/arg32;
  //F64vec8 ret64 = _mm512_cvtpslo_pd(arg32inv);
  //return ret64;
}

// construct weights for two particles based on their canonical positions
// 
// each particle has 8 weights, one for each corner of its cell
//
void construct_weights_for_2pcls(F64vec8 weights[2], F64vec8 X0)
{
  const F64vec8 X1 = F64vec8(1.) - X0;
  
  // construct weights for particles
  //
  F64vec8 wx[2], wy[2], wz[2];
  // need versions with broadcasted component,
  // but doing via swizzle works only for the x component
  //wsx[0] = ws[0].aaaa();
  //wsx[1] = ws[1].aaaa();
  const __mmask8 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
  const __mmask8 rmask_00110011 = _mm512_int2mask(0xcc); // mask=11001100
  const __mmask8 rmask_00001111 = _mm512_int2mask(0xf0); // mask=11110000
  // construct weights for first particle
  //
  // weights[0][0] = X0[0]*X0[1]*X0[2]; // weight000
  // weights[0][1] = X0[0]*X0[1]*X1[2]; // weight001
  // weights[0][2] = X0[0]*X1[1]*X0[2]; // weight010
  // weights[0][3] = X0[0]*X1[1]*X1[2]; // weight011
  // weights[0][4] = X1[0]*X0[1]*X0[2]; // weight100
  // weights[0][5] = X1[0]*X0[1]*X1[2]; // weight101
  // weights[0][6] = X1[0]*X1[1]*X0[2]; // weight110
  // weights[0][7] = X1[0]*X1[1]*X1[2]; // weight111
  //
  {
    // weights for first particle use low vectors of X0 and X1
    //wx[0] = cat_low_halves(wsx[0],wsx[1]);
    wx[0] = broadcast_mask_blend(X0[0],X1[0],rmask_00001111);
    wy[0] = broadcast_mask_blend(X0[1],X1[1],rmask_00110011);
    wz[0] = broadcast_mask_blend(X0[2],X1[2],rmask_01010101);
    weights[0] = wx[0]*wy[0]*wz[0];
  }
  // construct weights for second particle
  //
  // weights[1][0] = X0[4]*X0[5]*X0[6]; // weight000
  // weights[1][1] = X0[4]*X0[5]*X1[6]; // weight001
  // weights[1][2] = X0[4]*X1[5]*X0[6]; // weight010
  // weights[1][3] = X0[4]*X1[5]*X1[6]; // weight011
  // weights[1][4] = X1[4]*X0[5]*X0[6]; // weight100
  // weights[1][5] = X1[4]*X0[5]*X1[6]; // weight101
  // weights[1][6] = X1[4]*X1[5]*X0[6]; // weight110
  // weights[1][7] = X1[4]*X1[5]*X1[6]; // weight111
  //
  {
    // weights for second particle use high vectors of X0 and X1
    //wx[1] = cat_hgh_halves(wsx[0],wsx[1]);
    wx[1] = broadcast_mask_blend(X0[4],X1[4],rmask_00001111);
    wy[1] = broadcast_mask_blend(X0[5],X1[5],rmask_00110011);
    wz[1] = broadcast_mask_blend(X0[6],X1[6],rmask_01010101);
    weights[1] = wx[1]*wy[1]*wz[1];
  }
}

inline F64vec8 sample_field(
  F64vec8 weights,
  F64vec8 field_components[8])
{
  F64vec8 fields = F64vec8(0.);
  #pragma unroll
  for(int c=0; c<8; c++)
    fields += F64vec8(weights[c])*field_components[c];
  return fields;
}

inline F64vec8 compute_uvg_for_2pcls(F64vec8 uorig, F64vec8 B, F64vec8 E, F64vec8 qdto2mc)
{
  /* F64vec8::aaaa() etc is (defectively) not declared const */
  /*const*/ F64vec8 Om = qdto2mc*B;
  const F64vec8 Omx = Om.aaaa();
  const F64vec8 Omy = Om.bbbb();
  const F64vec8 Omz = Om.cccc();
  // if pushing more than 2 particles at a time, could make
  // use of reduction intrinsics (_mm512_mask_reduce_add_pd)
  // and could calculate the reciprocal more efficiently.
  const F64vec8 omsq_p1 = F64vec8(1.) + Omx*Omx + Omy*Omy + Omz*Omz;
  // is there a faster way to implement the reciprocal?
  const F64vec8 denom = reciprocal_F64vec8(omsq_p1);
  // solve the position equation
  /*const*/ F64vec8 ut = uorig + qdto2mc*E;
  const F64vec8 udotOm =
      ut.aaaa()*Omx
    + ut.bbbb()*Omy
    + ut.cccc()*Omz;

  // solve the velocity equation 
  //
  return (ut + cross_product(ut,Om)+udotOm*Om)*denom;
}

// test pushing all particles in a mesh cell
void test_push_pcls_in_cell_AoS_localized_vectorization()
{
  dfloat dto2 = usample();
  F64vec8 qdto2mc = F64vec8(usample());

  F64vec8* field_components = (F64vec8*) ::field_components;

  F64vec8 cellstart;
  F64vec8 dx_inv;
  for(int i=0;i<3;i++)
  {
    cellstart[i]=sample_openu();
    dx_inv[i]=sample_openu();
  }
  for(int i=0;i<3;i++)
  {
    cellstart[i+4]=cellstart[i];
    dx_inv[i+4]=dx_inv[i];
  }
  #pragma omp parallel num_threads(2)
  {
    // iterate over all particles in this mesh cell
    //
    const Time start = get_time();
    #pragma omp for
    for(int pi=0;pi<NUMPCLS;pi+=2)
    {
      SpeciesPcl* pcl = &pcls[pi];
      // because SpeciesPcl fits in cache line:
      ASSUME_ALIGNED(pcl);

      // gather position and velocity data from particles
      //
      F64vec8 pcl0 = *(F64vec8*)&pcl[0];
      F64vec8 pcl1 = *(F64vec8*)&pcl[1];
      const F64vec8 xorig = cat_low_halves(pcl0,pcl1);
      F64vec8 xavg = xorig;
      const F64vec8 uorig = cat_hgh_halves(pcl0,pcl1);

      // sample field for both particles
      //
      // compute canonical position.
      //
      // for each coordinate direction,
      // Xleft gives fraction of cell to the left of position and
      // Xrght gives fraction of cell to the right of position.
      const F64vec8 X = dx_inv*(xavg-cellstart);
      //
      // use canonical coordinates to interpolate field
      //
      F64vec8 weights[2];
      construct_weights_for_2pcls(weights, X);
      F64vec8 fields[2];
      // sample fields for first particle
      fields[0] = sample_field(weights[0],field_components);
      // sample fields for second particle
      fields[1] = sample_field(weights[1],field_components);
      const F64vec8 B = cat_low_halves(fields[0],fields[1]);
      const F64vec8 E = cat_hgh_halves(fields[0],fields[1]);

      // use sampled field to push particle block
      //
      const F64vec8 uavg = compute_uvg_for_2pcls(uorig, B, E, qdto2mc);
      // update average position
      xavg = xorig + uavg*dto2;

      // update state of particle (assuming this is last iteration)
      {
        const F64vec8 xnew = xavg+xavg - xorig;
        const F64vec8 unew = uavg+uavg - uorig;
        //
        // copy xnew and unew into particles
        //
        //copy012to012(pcl0,xnew);
        //copy456to012(pcl1,xnew);
        //copy012to456(pcl0,unew);
        //copy456to456(pcl1,unew);
        //
        const F64vec8 pcl0new = cat_low_halves(xnew, unew);
        const F64vec8 pcl1new = cat_hgh_halves(xnew, unew);
        copy012and456(pcl0,pcl0new);
        copy012and456(pcl1,pcl1new);
        //
        // could save using no-read stores( _mm512_storenr_pd),
        // but we just read this, so presumably it is still in cache.
        _mm512_storenr_pd(&pcl[0], pcl0);
        _mm512_storenr_pd(&pcl[1], pcl1);
      }
    }
    report_time(start);
  }
}

// This is declared noinline in order to force vectorization the way I want.
//
void sample_field(
  dfloat fields[2*D],
  dfloat weights[8],
  dfloat field_components[8][2*D])
{
  ASSUME_ALIGNED(fields);
  ASSUME_ALIGNED(weights);
  ASSUME_ALIGNED(field_components);
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
  printf("# intelligent would make no difference; for some reason\n");
  printf("# this is a bit slower but in single precision is faster.\n");
  printf("# \n");
  test_push_pcls_in_cell_SoA_vectorized();
  test_push_pcls_in_cell_SoA_vectorized();
  printf("#\n");
  printf("# the only real difference between the following and SoA\n");
  printf("# vectorization is that data must be gathered at the beginning\n");
  printf("# of the loop and scattered at the end, but unfortunately those\n");
  printf("# each seem to take as long as the push itself.\n");
  printf("#\n");
  test_push_pcls_in_cell_AoS_scatter_gather_vectorization();
  test_push_pcls_in_cell_AoS_scatter_gather_vectorization();
  printf("#\n");
  printf("# This is over 2 times slower than SoA            \n");
  printf("# vectorization, partly because only 75%          \n");
  printf("# vectorization can be expected for physical      \n");
  printf("# vectors and partly because dot product          \n");
  printf("# and reciprocal are computed with only 25%       \n");
  printf("# vectorization; I estimate that pushing 8        \n");
  printf("# particles at a time would allow this to be 1.75 \n");
  printf("# times slower than SoA vectorization.            \n");
  printf("#\n");
  test_push_pcls_in_cell_AoS_localized_vectorization();
  test_push_pcls_in_cell_AoS_localized_vectorization();
  printf("# \n");
  printf("# move particles stopping at boundary of mesh cell:\n");
  push_pcls_in_cell_SoA_stopping_at_face();
  push_pcls_in_cell_SoA_stopping_at_face();
  printf("#\n");
  printf("# move particles in 8-particle transposable blocks\n");
  printf("# stopping at boundary of mesh cell:\n");
  push_SoA_blocks_stopping_at_face();
  push_SoA_blocks_stopping_at_face();
  printf("#\n");
  printf("#\n");
  printf("# demonstrating that the time required to query the time is negligible:\n");
  printf("#\n");
  donothing_in_parallel();
  donothing_in_parallel();
  donothing_in_serial();
  donothing_in_serial();
}
