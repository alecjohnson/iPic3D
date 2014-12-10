#ifndef OutputWrapperTXT_h
#define OutputWrapperTXT_h
class Particles3Dcomm;
#include "Alloc.h"
#include "Setting.h"
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
    void write_velocity_distribution(int cycle, const Particles3Dcomm* pcls);
    void append_to_satellite_traces(const Grid3DCU *grid,
        const_arr3_double Bx, const_arr3_double By, const_arr3_double Bz,
        const_arr3_double Ex, const_arr3_double Ey, const_arr3_double Ez,
        const_arr4_double Jxs,const_arr4_double Jys,const_arr4_double Jzs,
        const_arr4_double rhons);
    void append_conserved_quantities(int cycle,
          double total_energy, double bogus_momentum,
          double Eenergy, double Benergy, double gas_energy);
  private:
    void append_to_velocity_distribution(int cycle,
      int is, double maxVel, long long *VelocityDist);
  private:
    // filenames for text-based output files
    string filename_cqsat;
    string filename_cq;
    string filename_ds;
    const Collective   *col;
    const VCtopology3D *vct;
    const Grid3DCU     *grid;
  private:
    void init_output_files();
  public:
    OutputWrapperTXT(const Setting& setting) :
      col(&setting.col()),
      vct(&setting.vct()),
      grid(&setting.grid())
    { init_output_files(); }
};

#endif // OutputWrapperTXT_h
