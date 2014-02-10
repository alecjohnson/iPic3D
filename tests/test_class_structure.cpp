#include "../include/debug.h"
#include "../include/Alloc.h"

#if 0 // pseudo-code for particle processing

 Instantiate the following:

   T   threads
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
  void CalculateMoments() = 0;
  void MovePcls(double dt) = 0;
};

class Grid
{
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
  array3<PclList> pclGrid;
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

template<class Pcl>
class PclList_AoS<Pcl>
{
  array1<Pcl> pclList;
}

// 512 bytes (cache-line-sized)
class SpeciesPcl
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
}

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
