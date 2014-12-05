
#include "MPIdata.h"
#include "iPic3D.h"
#include "debug.h"
#include "TimeTasks.h"
#include <stdio.h>

using namespace iPic3D;

int main(int argc, char **argv) {

 MPIdata::init(&argc, &argv);
 Parameters::init_parameters();
 {
  iPic3D::c_Solver solver(argc, argv);
  solver.Init();

  timeTasks.resetCycle();
  solver.computeMoments();
  for (int i = solver.FirstCycle(); i <= solver.FinalCycle(); i++)
  {
    if (solver.is_rank0())
      printf(" ======= Cycle %d ======= \n",i);

    timeTasks.resetCycle();
    solver.advanceEfield();
    solver.moveParticles();
    solver.advanceBfield();
    solver.computeMoments();
    solver.WriteOutput(i);
    // print out total time for all tasks
    timeTasks.print_cycle_times(i);
  }

  solver.Finalize();
 }
 // close MPI
 MPIdata::instance().finalize_mpi();

 return 0;
}
