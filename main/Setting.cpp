#include "mpi.h"
#include "Setting.h"
#include "Collective.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "MPIdata.h"

Setting::Setting(int argc, const char** argv)
{
  _col = new Collective(argc, argv);
  _vct = new VCtopology3D(*_col);
  // initialize the virtual cartesian topology 
  _vct->setup_vctopology(MPI_COMM_WORLD);
  _grid = new Grid3DCU(_col, _vct);

  // Print the initial settings to stdout and a file
  if(_vct->is_rank0())
  {
    MPIdata::instance().Print();
    _vct->Print();
    _col->Print();
    _col->save();
  }
}

Setting::~Setting()
{
  delete _grid;
  delete _vct;
  delete _col;
}
