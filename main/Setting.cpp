#include "Setting.h"

Setting::Setting(int argc, const char** argv)
{
  col = new Collective(argc, argv);
  vct = new VCtopology3D(*col);
  // initialize the virtual cartesian topology 
  vct->setup_vctopology(MPI_COMM_WORLD);
  grid = new Grid3DCU(col, vct);

  // Print the initial settings to stdout and a file
  if(vct->is_rank0())
  {
    MPIdata::instance().Print();
    vct->Print();
    col->Print();
    col->save();
  }
}

Setting::~Setting()
{
  delete grid;
  delete vct;
  delete col;
}
