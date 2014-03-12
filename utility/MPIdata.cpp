#ifndef NO_MPI
  #include <mpi.h>
#endif
#include <stdio.h>
#include <assert.h>
#include "MPIdata.h"

// code to check that init() is called before instance()
//
// no need for this to have more than file scope
int MPIdata::rank=-1;
int MPIdata::nprocs=-1;
static bool MPIdata_is_initialized=false;
bool MPIdata_assert_initialized()
{
  assert(MPIdata_is_initialized);
  return true;
}

MPIdata& MPIdata::instance()
{
  // This is executed on the first call to check that
  // MPIdata has first been initialized.
  static bool check = MPIdata_assert_initialized();
  static MPIdata* instance = new MPIdata;
  // After the first call, this is the only line
  // that is actually executed.
  return *instance;
}

void MPIdata::init(int *argc, char ***argv) {
  assert(!MPIdata_is_initialized);

 #ifdef NO_MPI
  rank = 0;
  nprocs = 1;
 #else // NO_MPI
  /* Initialize the MPI API */
  MPI_Init(argc, argv);

  /* Set rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Set total number of processors */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 #endif // NO_MPI

  MPIdata_is_initialized = true;
}

void MPIdata::finalize_mpi() {
 #ifndef NO_MPI
  MPI_Finalize();
 #endif
}

void MPIdata::Print(void) {
  printf("\nNumber of processes = %d\n"
         "-------------------------\n\n", get_nprocs());
}

// extern MPIdata *mpi; // instantiated in iPIC3D.cpp

