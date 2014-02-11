#include "debug.h"
#include "assert.h"
#include "ipicmath.h" // for pow2roundup
#include "ipicdefs.h" // for pfloat
#include "arraysfwd.h"
#define NO_NEW
#include "Alloc.h"

using namespace iPic3D;

// Prototype for Moment-Implicit Method (MIM),
// i.e. Implicit Moment Method (IMM).
// IMM is the traditional name; the name MIM
// was evidently later introduced, perhaps on the
// grounds that in the IMM discretization of the
// Maxwell-Kinetic system the one term that is
// *not* implicit is the "pressure" moments.
// A full name that really says what it is
// would be the '"Explicit Moment" Implicit Method'.
// MIM here is referenced as "Mi".

#if 0

   C+C communication buffers (6 <= C <= 26)
       designated as recv_buffer[direction] and send_buffer[direction].
       C=6 if communicating to face-neighbor processes and C=26
       if communicating diagonally as well.
   T*T thread-exchange buffers
       designated as exch_buffer[sending_thread][recving_thread].
   T*C thread-specific sending communication buffers
       designated as thread_buffer[sending_thread][direction].

   "Sorting buffer" refers to a communication or exchange buffer.

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

class Vector;
//
struct Scalar
{
  virtual void grad(Vector&)=0;
  virtual void laplacian(Vector&)=0;
};

struct Vector
{
  virtual Scalar& fetch_1()=0;
  virtual Scalar& fetch_2()=0;
  virtual Scalar& fetch_3()=0;
  virtual void div(Scalar&) = 0;
  virtual void laplacian(Vector&) = 0;
};

// SymTensor2 should derive from this.
//
struct Tensor2
{
  virtual Scalar& fetch_11()=0;
  virtual Scalar& fetch_12()=0;
  virtual Scalar& fetch_13()=0;
  virtual Scalar& fetch_21()=0;
  virtual Scalar& fetch_22()=0;
  virtual Scalar& fetch_23()=0;
  virtual Scalar& fetch_31()=0;
  virtual Scalar& fetch_32()=0;
  virtual Scalar& fetch_33()=0;
  virtual void div(Vector&)=0;
};

struct SymTensor2: public Tensor2
{
  virtual Scalar& fetch_11()=0;
  virtual Scalar& fetch_12()=0;
  virtual Scalar& fetch_13()=0;
  virtual Scalar& fetch_21(){return fetch_12();}
  virtual Scalar& fetch_22()=0;
  virtual Scalar& fetch_23()=0;
  virtual Scalar& fetch_31(){return fetch_31();}
  virtual Scalar& fetch_32(){return fetch_32();}
  virtual Scalar& fetch_33()=0;
  virtual void div(Vector&)=0;
};

// Abstract base class for moments of the MIM.
// Concrete representation is on a specific mesh.
struct MiMoments
{
  virtual Scalar& fetch_f()=0;
  virtual Vector& fetch_fu()=0;
  virtual SymTensor2& fetch_fuu()=0;
};

struct EMfield
{
  virtual Vector& fetch_B()=0;
  virtual Vector& fetch_E()=0;
};

// abstract base class that defines the field solver interface.
//
struct MiFieldSolver_interface
{
  virtual void AdvanceField(double dt) = 0;
};

// This field solver will assume periodic boundary conditions
//
class MiSpectralFieldSolver: public MiFieldSolver_interface
{
  void AdvanceField(double dt)
  {
    dprint(dt);
  }
};

// This field solver will handle general boundary conditions
// but will be only second-order accurate
//
struct MiFieldSolver: public MiFieldSolver_interface
{
  void AdvanceField(double dt)
  {
    dprint(dt);
  }
};

// skeleton. This needs to exchange
// fields and moments with MiFieldSolver_interface.
class MiKineticSolver_interface
{
 public:
  virtual void calculate_moments(MiMoments&)const = 0;
  // should be declared 
  virtual void move_pcls(const EMfield&, double dt) = 0;
};

// should support access of the form pclList[i].x,
// which could in fact also be supported
// with an underlying SoA representation.
//
template<class type>
class List
{
 private: // members
  type* list;
  int size; // number of particles in list
  int maxsize; // maximum number of particles
 public: // access
  void pop()
  {
    assert(size>0);
    size--;
  }
  void push_back(const type& element)
  {
    int newsize = size+1;
    resize_if_needed(newsize);
    list[size] = element;
    size=newsize;
  }
  inline type& operator[](int i)
  {
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_lt(i, size);
    #endif
    ALIGNED(list);
    return list[i];
  }
  void delete_element(int i)
  {
    // works even for last particle,
    // though pop would be faster in that case.
    // 
    // why doesn't this compile?
    //this->operator[i] = list[--size];
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_lt(i, size);
    #endif
    ALIGNED(list);
    list[i] = list[--size];
  }
 public: // memory
  ~List()
  {
    AlignedFree(list);
  }
  List():
    list(0),
    size(0),
    maxsize(0)
  {}
  // truncate and reallocate list
  void alloc(int in)
  {
    maxsize = in;
    free(list);
    list = AlignedAlloc(type,maxsize);
    size = 0;
  }
  void resize_if_needed(int needed_size)
  {
    if(needed_size>maxsize)
    {
      int newmaxsize = pow2roundup(needed_size);
      resize(newmaxsize);
    }
  }
  // change maximum size without deleting elements
  void resize(int newmaxsize)
  {
    type* oldList = list;
    assert(size < newmaxsize);
    maxsize = newmaxsize;
    for(int i=0;i<size;i++)
    {
      list[i] = oldList[i];
    }
    AlignedFree(oldList);
  }
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
  PclList& fetch(int i)
  {
    return pclGrid.fetch(i);
  }
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
  // eliminate (temporary convenience for testing)
  PclGrid(int nxc_, int nyc_, int nzc_)
   :Grid(nxc_, nyc_, nzc_, 1, 1, 1, 0, 0, 0),
    pclGrid(nxc,nyc,nzc)
  {}
  // resize all particle lists to hold as many as npc particles
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

struct PclRef
{
  // grid element of list where particle resides
  // (assumes PclGrid has underlying 1D array)
  int cell;
  // index of particle in list
  int pcl;
};

class PclRefList: public List<PclRef>
{
};

template<class PclList>
class PclSolver:
  public PclGrid<PclList>,
  public MiKineticSolver_interface
{
  // particle sorting buffers
  //
  // number of directions to communicate
  int C; // perhaps 6 or 26
  // number of threads per process
  int T;
  array1<PclList> send_buffers; // C
  array1<PclList> recv_buffers; // C
  array2<PclRefList> tsend_buffers; // TxC
  array2<PclRefList> texch_buffers; // TxT
 public:
  PclSolver(int C_, int T_)
   :C(C_),
    T(T_),
    send_buffers(C),
    recv_buffers(C),
    tsend_buffers(T,C),
    texch_buffers(T,T)
  {
  }
  void calculate_moments(MiMoments&)const
  {
    dprint("gothere");
  }
  void move_pcls(const EMfield&, double dt)
  {
    dprint(dt);
  }
  void sort_particles()
  {
    
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
class MiPcl
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
class MiPcl3
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
//  //PclGrid<List<SpeciesPcl>> pclGrid;
//  List<SpeciesPcl> pclList;
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
  PclGrid<List<SpeciesPcl> > pclGrid(2,3,4);
  pclGrid.alloc(2); // at most 2 particles per cell
  SpeciesPcl pcl;
  pcl.x=1.234;
  pclGrid.fetch(1,2,3).push_back(pcl);
  pcl.x=2.345;
  pclGrid.fetch(1,2,3).push_back(pcl);
  List<SpeciesPcl>& pclList = pclGrid.fetch(1,2,3);
  //SpeciesPcl& pcl0 = pclList[0];
  printf(
    "pclGrid[%d][%d][%d][%d].x = %g\n",
    1,2,3,0,pclList[0]);
  printf(
    "pclGrid[%d][%d][%d][%d].x = %g\n",
    1,2,3,1,pclList[1]);
}

void test_sort()
{
}

int main()
{
  printf("hello world\n");
  test_PclGrid();
}
