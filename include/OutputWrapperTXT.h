#ifndef OutputWrapperTXT_h
#define OutputWrapperTXT_h
#include "Alloc.h"
#include<string>

using std::string;
// ===
// OutputWrapperFPP: output wrapper for text file output
//
//   This class should provide a mechanism to avoid 
//   repeated opening and closing of the same file.
// ===

class OutputWrapperTXT
{
  public:
    static bool do_write_virtual_satellite_traces();
    static bool do_write_velocity_distribution();
    static int get_number_of_distribution_bins(){return 1000;}
  private:
    // filenames for text-based output files
    string filename_cqsat;
    string filename_cq;
    string filename_ds;
  private:
    void init_output_files(
      Collective    *col,
      VCtopology3D  *vct);
  public:
    OutputWrapperTXT(
      Collective    *col,
      VCtopology3D  *vct)
    { init_output_files(col,vct); }
};

#endif // OutputWrapperTXT_h
