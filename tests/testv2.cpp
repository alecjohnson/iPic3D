#include <omp.h>
#include "stdio.h"
#include <stdlib.h> // rand()
#include <sys/time.h>
#include <limits.h> // RAND_MAX
#include <string.h> // memcpy
#include "../utility/debug.cpp"

#define ALIGNMENT 64
#define ALLOC_ALIGNED __attribute__((aligned(ALIGNMENT)));
#define ASSERT_ALIGNED(X) __assume_aligned(X, ALIGNMENT)

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

inline double time_msec()
{
  // this seems more reliable; gettimeofday even results in negative times
  return 1000.*omp_get_wtime();
  static struct timeval tv;

  gettimeofday(&tv, NULL);

  return (tv.tv_sec + tv.tv_usec * (double) 1e-3);
}

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

double usample()
{
  const double RAND_MAX_inv = 1./RAND_MAX;
  return rand()*RAND_MAX_inv;
}

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
    for(int i=0;i<3;i++)
    {
      x[i] = usample();
      u[i] = usample();
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

void copy_SoA_to_AoS_particles()
{
  #pragma omp parallel
  {
  const double start = time_msec();
  #pragma omp for simd
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
  }
}

void copy_AoS_to_SoA_particles()
{
  #pragma omp parallel
  {
  const double start = time_msec();
  #pragma omp for simd
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
  }
}

void copy_AoS_particles_to_another_array()
{
  #pragma omp parallel
  {
  const double start = time_msec();
  #pragma omp for
  for(int p=0;p<NUMPCLS;p++)
  {
    ASSERT_ALIGNED(&pcls2[0]);
    ASSERT_ALIGNED(&pcls[0]);
    pcls2[p] = pcls[p];
  }
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
  }
}

void copy_AoS_particles_to_another_array_via_memcpy()
{
  const double start = time_msec();
  memcpy(&pcls2[0],&pcls[0],sizeof(pcls[0])*NUMPCLS);
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
}

void copy_AoS_to_AoS_particle_blocks()
{
  #pragma omp parallel
  {
  const double start = time_msec();
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
  }
}

void copy_AoS_to_AoS_particle_blocks_via_memcpy()
{
  // How to parallelize memcpy?
  const double start = time_msec();
  memcpy(&pclBlocks[0],&pcls[0],sizeof(pcls[0])*NUMPCLS);
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
}

void initialize_data()
{
  for(int c=0;c<8;c++)
  for(int i=0;i<2*D;i++)
    field_components[c][i] = usample();

  // initialize particles
  //
  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].init_random(i);
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
  const double start = time_msec();
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
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
  const double start = time_msec();
  #pragma omp for simd
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
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
        const double start = time_msec();
        #pragma omp for simd
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
            // fraction of the distance from the right of the cell
            const double w1x = dx_inv*(xavg - cellstartx);
            const double w1y = dy_inv*(yavg - cellstarty);
            const double w1z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the left
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
        const double end = time_msec();
        printf("(in thread %d) %s took %g ms.\n",
          omp_get_thread_num(), __func__, end - start);
    }
}

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
      const double start = time_msec();
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
            // fraction of the distance from the right of the cell
            const double w1x = dx_inv*(xavg - cellstartx);
            const double w1y = dy_inv*(yavg - cellstarty);
            const double w1z = dz_inv*(zavg - cellstartz);
            // fraction of distance from the left
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
      const double end = time_msec();
      printf("(in thread %d) %s took %g ms.\n",
        omp_get_thread_num(), __func__, end - start);
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
  const double start = time_msec();
  //simd accelerates execution by a factor of 3
  #pragma omp for simd
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
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
  const double start = time_msec();
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
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
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
  const double start = time_msec();
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
}

void donothing_in_parallel()
{
  #pragma omp parallel
  {
  const double start = time_msec();
  const double end = time_msec();
  printf("(in thread %d) %s took %g ms.\n",
    omp_get_thread_num(), __func__, end - start);
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
  printf("#\n");
  printf("#\n");
  printf("# demonstrating that the time required to query the time is negligible:\n");
  printf("#\n");
  donothing_in_parallel();
  donothing_in_parallel();
  donothing_in_serial();
  donothing_in_serial();
}
