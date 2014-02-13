#include "stdio.h"
#include <stdlib.h> // rand()
#include <limits.h> // RAND_MAX

#define ALIGNMENT 64
#define ALLOC_ALIGNED __attribute__((aligned(ALIGNMENT)));
#define ASSERT_ALIGNED(X) __assume_aligned(X, ALIGNMENT)

// const int D=4;
// const int NPB=2; // number of particles in a block
#define D (4)
#define NPB (2)
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
  double field_components[8][2*D])__attribute__((noinline));

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
SpeciesPcl pcls[NUMPCLS];

void initialize_data()
{
  for(int c=0;c<8;c++)
  for(int i=0;i<2*D;i++)
    field_components[c][i] = usample();

  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].init_random(i);
  }
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
    #pragma simd collapse(2)
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
    double uavg[D];
    {
      double Om[D];
      double denom;
      for(int i=0;i<D;i++)
      {
        Om[i] = qdto2mc*B[i];
      }
      {
        double omsq_p1 = 1. + Om[0]*Om[0] + Om[1]*Om[1] + Om[2]*Om[2];
        denom = 1.0f/float(omsq_p1);
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
      #pragma simd collapse(2)
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
}

void test_push_pcls_in_cell_SoA()
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
      #pragma unroll
      for(int c=0; c<8; c++)
      #pragma unroll
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
}

// test pushing all particles in a mesh cell
void test_push_pcls_in_cell_vectorized_AoS()
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
    #pragma simd
    for(int i=0;i<D;i++)
    {
      xorig[p][i] = pcl[p].get_x(i);
      uorig[p][i] = pcl[p].get_u(i);
    }
    #pragma simd collapse(2)
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
    // if(vectorized_w)
    #pragma simd collapse(2)
    for(int p=0; p<NPB; p++)
    for(int i=0;i<D;i++)
    {
      ws[1][p][i] = dx_inv[i]*(xavg[p][i]-cellstart[i]);
      ws[0][p][i] = 1.-ws[1][p][i];
    }
    for(int p=0; p<NPB; p++)
    {
      double fields[NPB][D]ALLOC_ALIGNED;
      double weights[8]ALLOC_ALIGNED;
      //if(vectorized_w)
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
      //else
      {
        double w[2][D]ALLOC_ALIGNED;
        #pragma simd
        for(int i=0;i<D;i++)
        {
          w[1][i] = dx_inv[i]*(xavg[p][i]-cellstart[i]);
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
      sample_field(fields[p],weights,field_components);
      #pragma simd
      // scatter field data for vectorized push
      for(int i=0;i<D;i++)
      {
        B[p][i] = fields[p][i];
        E[p][i] = fields[p][D+i];
      }
    }

    // use sampled field to push particle block
    //
    double uavg[NPB][D]ALLOC_ALIGNED;
    {
      double Om[NPB][D]ALLOC_ALIGNED;
      double denom[NPB]ALLOC_ALIGNED;
      #pragma simd
      for(int i=0;i<NPB*D;i++)
      {
        Om[0][i] = qdto2mc*B[0][i];
      }
      for(int p=0;p<NPB;p++)
      {
        double omsq_p1 = 1.
                     + Om[p][0] * Om[p][0]
                     + Om[p][1] * Om[p][1]
                     + Om[p][2] * Om[p][2];
        denom[p] = 1.0f/omsq_p1;
      }
      double ut[NPB][D]ALLOC_ALIGNED;
      double udotOm[NPB]ALLOC_ALIGNED;
      // solve the position equation
      //#pragma simd collapse(2)
      //for(int p=0;p<NPB;p++)
      //for(int i=0;i<D;i++)
      //{
      //  ut[p][i] = uorig[p][i] + qdto2mc*E[p][i];
      //}
      #pragma simd
      for(int i=0;i<NPB*D;i++)
      {
        ut[0][i] = uorig[0][i] + qdto2mc*E[0][i];
      }
      for(int p=0;p<NPB;p++)
      {
        udotOm[p] = ut[p][0]*Om[p][0] + ut[p][1]*Om[p][1] + ut[p][2]*Om[p][2];
      }
      // solve the velocity equation 
      //
      #pragma simd // how do I tell it to recognize the swizzle?
      for(int p=0;p<NPB;p++)
      {
        uavg[p][0] = (ut[p][0] + (ut[p][1] * Om[p][2] - ut[p][2] * Om[p][1] + udotOm[p] * Om[p][0])) * denom[p];
        uavg[p][1] = (ut[p][1] + (ut[p][2] * Om[p][0] - ut[p][0] * Om[p][2] + udotOm[p] * Om[p][1])) * denom[p];
        uavg[p][2] = (ut[p][2] + (ut[p][0] * Om[p][1] - ut[p][1] * Om[p][0] + udotOm[p] * Om[p][2])) * denom[p];
        if(D==4) // to help with recognizing swizzle
        uavg[p][3] = (ut[p][3] + (ut[p][3] * Om[p][3] - ut[p][3] * Om[p][3] + udotOm[p] * Om[p][3])) * denom[p];
      }
      // update average position
      #pragma simd collapse(2)
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

int main()
{
  test_alex();
  test_alex2();
  test_AoS();
  initialize_data();
  test_push_pcls_in_cell();
  test_push_pcls_in_cell_vectorized_AoS();
  test_push_pcls_in_cell_SoA();
}
