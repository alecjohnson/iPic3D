
#include <mpi.h>
#include <iomanip>
#include "iPic3D.h"
#include "debug.h"

using namespace iPic3D;

int main(int argc, char **argv) {

  iPic3D::c_Solver KCode;
  bool b_err = false;

  MPIdata::init(&argc, &argv);
  dprintf("gothere");
  KCode.Init(argc, argv);
  dprintf("gothere");

  for (int i = KCode.FirstCycle(); i < KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    if (!b_err) {
  dprintf("gothere");
      KCode.CalculateField();
  dprintf("gothere");
      b_err = KCode.ParticlesMover();
  dprintf("gothere");
    }

    if (b_err) {
      i = KCode.LastCycle() + 1;
    }

    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

  }

  KCode.Finalize();

  return 0;
}
