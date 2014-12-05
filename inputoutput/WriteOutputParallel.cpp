
#include <mpi.h>
#include "WriteOutputParallel.h"
#include "phdf5.h"
#include "CollectiveIO.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "VCtopology3D.h"
#include "errors.h"
#include <string>
#include <sstream>
#include <iomanip>
using std::string;

void WriteOutputParallel(
  CollectiveIO *col,
  VCtopology3D *vct,
  Grid3DCU *grid,
  Particles3Dcomm *part,
  Pmoments *pMoments,
  EMfields3D *EMf,
  int cycle)
{
#ifdef PHDF5
  string       grpname;
  string       dtaname;

  stringstream filenmbr;
  string       filename;

  /* ------------------- */
  /* Setup the file name */
  /* ------------------- */

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  /* ---------------------------------------------------------------------------- */
  /* Define the number of cells in the globa and local mesh and set the mesh size */
  /* ---------------------------------------------------------------------------- */

  int nxc_r = grid->get_nxc_r();
  int nyc_r = grid->get_nyc_r();
  int nzc_r = grid->get_nzc_r();

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { nxc_r,            nyc_r,            nzc_r };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  /* --------------------------------------- */
  /* Declare and open the parallel HDF5 file */
  /* --------------------------------------- */

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getComm());

  outputfile.CreatePHDF5file(L, dglob, dlocl, false);

  /* ------------------------ */
  /* Write the Electric field */
  /* ------------------------ */

  outputfile.WritePHDF5dataset("Fields", "Ex", EMf->ret_Exc(), nxc_r, nyc_r, nzc_r);
  outputfile.WritePHDF5dataset("Fields", "Ey", EMf->ret_Eyc(), nxc_r, nyc_r, nzc_r);
  outputfile.WritePHDF5dataset("Fields", "Ez", EMf->ret_Ezc(), nxc_r, nyc_r, nzc_r);

  /* ------------------------ */
  /* Write the Magnetic field */
  /* ------------------------ */

  outputfile.WritePHDF5dataset("Fields", "Bx", EMf->ret_Bxc(), nxc_r, nyc_r, nzc_r);
  outputfile.WritePHDF5dataset("Fields", "By", EMf->ret_Byc(), nxc_r, nyc_r, nzc_r);
  outputfile.WritePHDF5dataset("Fields", "Bz", EMf->ret_Bzc(), nxc_r, nyc_r, nzc_r);

  /* ----------------------------------------------- */
  /* Write the moments for each species */
  /* ----------------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream snmbr;
    snmbr << is;
    const string num = snmbr.str();

    // Charge Density
    outputfile.WritePHDF5dataset("Fields", string("Rho_")+num , pMoments->ret_rhocs(is), nxc_r, nyc_r, nzc_r);
    // Current
    outputfile.WritePHDF5dataset("Fields", string("Jx_")+num, pMoments->ret_Jxsc(is), nxc_r, nyc_r, nzc_r);
    outputfile.WritePHDF5dataset("Fields", string("Jy_")+num, pMoments->ret_Jysc(is), nxc_r, nyc_r, nzc_r);
    outputfile.WritePHDF5dataset("Fields", string("Jz_")+num, pMoments->ret_Jzsc(is), nxc_r, nyc_r, nzc_r);
  }

  outputfile.ClosePHDF5file();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}
