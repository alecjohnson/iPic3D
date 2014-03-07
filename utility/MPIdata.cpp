#include <mpi.h>
#include <assert.h>
#include "MPIdata.h"
#include "ompdefs.h" // for omp_get_max_threads

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

  /* Initialize the MPI API */
  MPI_Init(argc, argv);

  /* Set rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Set total number of processors */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MPIdata_is_initialized = true;
}

void MPIdata::finalize_mpi() {
  MPI_Finalize();
}

void MPIdata::Print(void) {
  printf("\n"
    "Number of processes = %d\n"
    "-------------------------\n"
    "Number of threads = %d\n"
    "-------------------------\n",
     get_nprocs(),
     omp_get_max_threads());
}

// extern MPIdata *mpi; // instantiated in iPIC3D.cpp

