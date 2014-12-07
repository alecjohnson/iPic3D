
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
  solver.run();
 }
 // close MPI
 MPIdata::instance().finalize_mpi();

 return 0;
}
