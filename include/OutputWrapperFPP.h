#ifndef OutputWrapperFPP_h
#define OutputWrapperFPP_h
// ===
// OutputWrapperFPP: output wrapper for file-per-process output
//
//   This class should provide a mechanism to avoid 
//   repeated opening and closing of the same file.
// ===
#include "ipic_fwd.h"
#include "PSKOutput.h"
#include "PSKhdf5adaptor.h"

using namespace PSK;

class OutputWrapperFPP
{
 private:
  const Setting       &setting;
  const Collective    *col;
  const VCtopology3D  *vct;
  const Grid3DCU      *grid;
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
  OutputWrapperFPP(const Setting& setting_) :
    setting(setting_),
    col(&setting.col()),
    vct(&setting.vct()),
    grid(&setting.grid())
  {}
  void init_output_files(
    const EMfields3D    *EMf,
    const Particles3D   *part);
  void append_output(const char* tag, int cycle);
  void append_output(const char* tag, int cycle, int sample);
};

#endif // OutputWrapperFPP_h
