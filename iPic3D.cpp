
#include <mpi.h>
#include <iomanip>
#include "iPic3D.h"
#include "debug.h"
#include "TimeTasks.h"

using namespace iPic3D;

int main(int argc, char **argv) {

 MPIdata::init(&argc, &argv);
 {
  iPic3D::c_Solver KCode;
  KCode.Init(argc, argv);

  timeTasks.resetCycle();
  KCode.CalculateMoments();
  for (int i = KCode.FirstCycle(); i < KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    timeTasks.resetCycle();
    KCode.CalculateField();
    KCode.ParticlesMover();
    KCode.CalculateB();
    KCode.CalculateMoments();
    KCode.WriteOutput(i);
    // print out total time for all tasks
    timeTasks.print_cycle_times(i);
  }

  KCode.Finalize();
 }
 // close MPI
 MPIdata::instance().finalize_mpi();

 return 0;
}
