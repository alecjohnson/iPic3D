#ifndef OutputWrapperFPP_h
#define OutputWrapperFPP_h
// ===
// OutputWrapperFPP: output wrapper for file-per-process output
//
//   This class should provide a mechanism to avoid having
//   repeatedly opening and closing the same file.
// ===
#include "ipicfwd.h"
#include "PSKOutput.h"
#include "PSKhdf5adaptor.h"

using namespace PSK;

class OutputWrapperFPP
{
 private:
  #ifndef NO_HDF5
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output
  #endif // NO_HDF5
  int cartesian_rank;
  string SaveDirName;
  string RestartDirName;
  string output_file;
 public:
  void init_output_files(
    Collective    *col,
    VCtopology3D  *vct,
    Grid3DCU      *grid,
    EMfields3D    *EMf,
    Particles3D   *part,
    int ns);
  void append_output(const char* tag, int cycle);
  void append_output(const char* tag, int cycle, int sample);
};

#endif // OutputWrapperFPP_h
