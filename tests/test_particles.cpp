#include "../include/debug.h"
#include "../include/assert.h"
#include "../include/ipicdefs.h" // for pfloat
#include "../include/arraysfwd.h"
#define NO_NEW
#include "../include/Alloc.h"

using namespace iPic3D;

#if 0

   C+C communication buffers (6 <= C <= 26)
       designated as recv_buffer[direction] and send_buffer[direction].
       C=6 if communicating to face-neighbor processes and C=26
       if communicating diagonally as well.
   T*T thread-exchange buffers
       designated as exch_buffer[sending_thread][recving_thread].
   T*C thread-specific sending communication buffers
       designated as thread_buffer[sending_thread][direction].

   "Sorting buffer" refers to a communication or exchange buffers.

 Move particles:

   #pragma omp parallel
   foreach iteration of implicit pusher:
     if sorting_with_each_time_step:
       assert: particles are already sorted by xavg
     foreach cell in subdomain of this thread:
       foreach particle in cell:
         a. interpolate field at xavg
         b. iterate particle advance and update xavg
         c. if final_iteration:
              set x = 2*xavg-x
              set u = 2*uavg-u
              // maintained for sorting and initial push:
              set xavg = x
         d. if sorting before next push and xavg is outside cell:
              if xavg is inside thread subdomain:
                if sorting_and_moving_simultaneously:
                  push particle to destination cell
                else:
                  append particle index to list of particles
                  ...to be gathered by recipient cell
              else if xavg is inside process subdomain:
                append particle index to list of particles to be
                ...pulled from this thread by recipient thread
              else: // xavg is in another process subdomain
                append particle index to list of particles to be
                ...pulled from this thread by recipient communication buffer
     if sorting_with_each_time_step:
       sort particles by xavg

 Sort as follows:

   #pragma omp (parallel) for
   1. Each thread identifies particles exiting its subdomain
      and places [index of] each exiting particle in an
      appropriate sorting buffer [already done in mover].

   done_communicating = false;
   while sorting buffers are not empty:
   {
     if communication send buffers are empty:
       done_communicating = true;

     if not done_communicating
     {
       #pragma omp barrier
       #pragma omp single
       2. Allocate communication buffers:
          a. Calculate needed size of all C buffers.
          b. Communicate buffer sizes.
          c. Resize receiving buffers if needed.
          d. Call MPI_Irecv on each recv_buffer.
          e. Resize each send_buffer if needed.
       #pragma omp (parallel) for
       3. Concatenate send buffers:
            Each thread copies the outgoing particles referenced in
            its private communication send buffers to a calculated
            segment in the corresponding common communication buffer.
       #pragma omp single
       4. Send communicated particles:
          Call MPI_Isend on each send_buffer.
     }
     #pragma omp (parallel) for
     5. Sort within each subdomain:
          if first_iteration:
            Each thread sorts the particles
              that move within its subdomain.
          else:
            Each thread walks through the particles in the exchange
            buffers it reads and transfers each particle to its
            correct mesh cell.
     if not done_communicating:
     {
       MPI_Waitall(...,recv_requests);
       6. Move particles from receive buffers to sorting buffers:
            (#pragma omp parallel)
            foreach recv_buffer:
              #pragma omp for
              foreach particle in recv_buffer:
                put particle in appropriate sorting buffer
     }
   }
#endif

// Proposed class structure for iPic3D.
//
// The purpose of this program is to test
// that the template structure and abstract base
// classes will compile and run properly.

// abstract base class that defines the field solver interface.
//
struct ImFieldSolver_interface
{
  virtual void AdvanceField(double dt) = 0;
};

// This field solver will assume periodic boundary conditions
//
class ImSpectralFieldSolver: public ImFieldSolver_interface
{
  void AdvanceField(double dt)
  {
    dprint(dt);
  }
};

// This field solver will handle general boundary conditions
// but will be only second-order accurate
//
struct ImFieldSolver: public ImFieldSolver_interface
{
  void AdvanceField(double dt)
  {
    dprint(dt);
  }
};

class ImKineticSolver_interface
{
 public:
  virtual void CalculateMoments() = 0;
  virtual void MovePcls(double dt) = 0;
};

class Grid
{
 protected:
  int nxc;
  int nyc;
  int nzc;
  double dx;
  double dy;
  double dz;
  double xStart;
  double yStart;
  double zStart;
 public:
  Grid(int nxc_, int nyc_, int nzc_,
    double dx_, double dy_, double dz_,
    double xStart_, double yStart_, double zStart_)
   :nxc(nxc_), nyc(nyc_), nzc(nzc_),
    dx(dx_), dy(dy_), dz(dz_),
    xStart(xStart_), yStart(yStart_), zStart(zStart_)
  {}
  void get_cell_for_position(
    int& xc, int& yc, int& zc,
    double x, double y, double z);
};

class FieldGrid: public Grid
{
 public:
  void divergence(
    array_ref3<double> out,
    const_array_ref3<double> xin,
    const_array_ref3<double> yin,
    const_array_ref3<double> zin);
};

template<class PclList>
class PclGrid: public Grid
{
 public: // will be private
  array3<PclList> pclGrid;
 public:
  PclList& fetch(int i, int j, int k)
  {
    return pclGrid[i][j][k];
  }
  PclGrid(int nxc_, int nyc_, int nzc_,
    double dx_, double dy_, double dz_,
    double xStart_, double yStart_, double zStart_)
   :Grid(nxc_, nyc_, nzc_, dx_, dy_, dz_,
      xStart_, yStart_, zStart_),
    pclGrid(nxc,nyc,nzc)
  {}
  PclGrid(int nxc_, int nyc_, int nzc_)
   :Grid(nxc_, nyc_, nzc_, 1, 1, 1, 0, 0, 0),
    pclGrid(nxc,nyc,nzc)
  {}
  // resize all particle lists to hold as much as npc particles
  void resize(int npc)
  {
    for(int i=0;i<nxc;i++)
    for(int j=0;j<nyc;j++)
    for(int k=0;k<nzc;k++)
      pclGrid[i][j][k].resize(npc);
  }
  void alloc(int npc)
  {
    for(int i=0;i<nxc;i++)
    for(int j=0;j<nyc;j++)
    for(int k=0;k<nzc;k++)
      pclGrid[i][j][k].alloc(npc);
  }
};

template<class PclList>
class PclSolver:
  public PclGrid<PclList>,
  public ImKineticSolver_interface
{
 public:
  void CalculateMoments(){}
  void MovePcls(double dt)
  {
    dprint(dt);
  }
};

// should support access of the form pclList[i].x,
// which could in fact also be supported
// with an underlying SoA representation.
//
template<class Pcl>
class PclList_AoS
{
 private: // members
  Pcl* pclList;
  int last; // index of last particle in list
  int maxsize; // maximum number of particles
 public: // access
  void pop()
  {
    assert(last>=0);
    last--;
    // return pclList[last--];
  }
  void push_back(const Pcl& pcl)
  {
    last++;
    pclList[last] = pcl;
    assert(last<maxsize);
  }
  inline Pcl& operator[](int i)
  {
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_le(i, last);
    #endif
    ALIGNED(pclList);
    return pclList[i];
  }
  void delete_pcl(int i)
  {
    // works even for last particle,
    // though pop would be faster in that case.
    // 
    // why doesn't this compile?
    //this->operator[i] = pclList[last--];
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_le(i, last);
    #endif
    ALIGNED(pclList);
    pclList[i] = pclList[last--];
  }
 public: // memory
  ~PclList_AoS()
  {
    AlignedFree(pclList);
  }
  PclList_AoS():
    pclList(0),
    last(-1),
    maxsize(0)
  {}
  // truncate and reallocate list
  void alloc(int in)
  {
    maxsize = in;
    free(pclList);
    pclList = AlignedAlloc(Pcl,maxsize);
    last = -1;
  }
  // change maximum size without deleting pcls
  void resize(int newmaxsize)
  {
    Pcl* oldPclList = pclList;
    assert(last+1 < newmaxsize);
    maxsize = newmaxsize;
    for(int i=0;i<=last;i++)
    {
      pclList[i] = oldPclList[i];
    }
    AlignedFree(oldPclList);
  }
};

// 512 bytes (cache-line-sized)
struct SpeciesPcl
{
  double x;
  double y;
  double z;
  double q;
  double u;
  double v;
  double w;
  long long ID;
};

// 512 bytes (cache-line-sized)
class ImPcl
{
  int cx;
  int cy;
  int cz;
  float dx;
  float dy;
  float dz;
  float xavg;
  float yavg;
  float zavg;
  float u;
  float v;
  float w;
  float m; // mass
  float qom; // charge-to-mass ratio
  long long ID;
};

// 1024 bits, for 3rd-order mover
//
class ImPcl3
{
  double x;
  double y;
  double z;
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  double u;
  double v;
  double w;
  double m; // mass
  double qom; // charge-to-mass ratio
  long long ID;
};

//int test_arrays()
//{
//  array1_int arr(5);
//  //int* arr = (int*)malloc(sizeof(int)*5);
//  for(int i=0;i<5;i++)
//    arr[i] = i;
//  for(int i=0;i<5;i++)
//    printf("arr[%d] = %d\n", i, arr[i]);
//  array3<int> arr3(4,3,2);
//  for(int i=0;i<4;i++)
//  for(int j=0;j<3;j++)
//  for(int k=0;k<2;k++)
//    arr3[i][j][k] = i+j+k;
//  for(int i=0;i<4;i++)
//  for(int j=0;j<3;j++)
//  for(int k=0;k<2;k++)
//    printf("arr3[%d][%d][%d] = %d\n",i,j,k, arr3[i][j][k]);
//}
//
//int test_PclList()
//{
//  // create grid of particles
//  //PclGrid<PclList_AoS<SpeciesPcl>> pclGrid;
//  PclList_AoS<SpeciesPcl> pclList;
//  pclList.alloc(5);
//  SpeciesPcl pcl0;
//  SpeciesPcl pcl1;
//  pclList.push_back(pcl0);
//  pclList.push_back(pcl1);
//  pclList[1].x = 1.234;
//  printf("pclList[%d].x=%g\n", 1, pclList[1].x);
//  pcl0.y = 2.345;
//  pclList[0] = pcl0;
//  printf("pclList[%d].y=%g\n", 0, pclList[0].y);
//
//  array1<SpeciesPcl> pcls(5);
//  pcls[1].x = 5.342;
//  printf("pcls[%d].x=%g\n", 1, pcls[1].x);
//}

int test_PclGrid()
{
  // create grid of particles
  PclGrid<PclList_AoS<SpeciesPcl> > pclGrid(2,3,4);
  pclGrid.alloc(2); // at most 2 particles per cell
  SpeciesPcl pcl;
  pcl.x=1.234;
  pclGrid.fetch(1,2,3).push_back(pcl);
  pcl.x=2.345;
  pclGrid.fetch(1,2,3).push_back(pcl);
  PclList_AoS<SpeciesPcl>& pclList = pclGrid.fetch(1,2,3);
  //SpeciesPcl& pcl0 = pclList[0];
  printf(
    "pclGrid[%d][%d][%d][%d].x = %g\n",
    1,2,3,0,pclList[0]);
  printf(
    "pclGrid[%d][%d][%d][%d].x = %g\n",
    1,2,3,1,pclList[1]);
}

int main()
{
  printf("hello world\n");
  test_PclGrid();
}
