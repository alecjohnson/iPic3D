#include "OutputWrapperTXT.h"
#include <fstream>

bool OutputWrapperTXT::do_write_virtual_satellite_traces(){return false;}
bool OutputWrapperTXT::do_write_velocity_distribution(){return false;}
int OutputWrapperTXT::get_number_of_distribution_bins(){return 1000;}
static int get_number_of_satellites_per_direction(){return 3;}

// implementation of a C++11 method
static string to_string(int n)
{
  stringstream n_ss;
  n_ss << n;
  return n_ss.str();
}

// initialize output files
OutputWrapperTXT::init_output_files(
  Collective    *col,
  VCtopology3D  *vct);
{
  filename_cq = col->getSaveDirName() + "/ConservedQuantities.txt";
  if (myrank == 0) {
    ofstream my_file(filename_cq.c_str());
    my_file.close();
  }
  // Distribution functions
  if(do_write_distribution_functions())
  {
    const int nDistributionBins = get_number_of_distribution_bins();
    filename_ds = col->getSaveDirName() + "/DistributionFunctions.txt";
    if (myrank == 0) {
      ofstream my_file(filename_ds.c_str());
      my_file.close();
    }
  }

  // set virtual satellite traces
  //
  if(do_write_virtual_satellite_traces())
  {
    filename_cqsat = col->getSaveDirName()
      + "/VirtualSatelliteTraces"
      + to_string(vct->get_rank())
      + ".txt";
    // if(myrank==0)
    ofstream my_file(filename_cqsat.c_str(), fstream::binary);
    const int nsat = get_number_of_satellites_per_direction();
    const int nx0 = grid->get_nxc_r();
    const int ny0 = grid->get_nyc_r();
    const int nz0 = grid->get_nzc_r();
    for (int isat = 0; isat < nsat; isat++)
    for (int jsat = 0; jsat < nsat; jsat++)
    for (int ksat = 0; ksat < nsat; ksat++)
    {
      int i1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
      int i2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
      int i3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
      my_file << grid->getXC(i1) << "\t"
              << grid->getYC(i2) << "\t"
              << grid->getZC(i3) << endl;
    }
    my_file.close();
  }
}

void OutputWrapperTXT::append_conserved_quantities(
      double total_energy, double bogus_momentum,
      double Eenergy, double Benergy, double gas_energy)
{
  if (myrank == 0) {
    ofstream my_file(filename_cq.c_str(), fstream::app);
    my_file << cycle << "\t"
      << "\t" << total_energy
      << "\t" << bogus_momentum
      << "\t" << Eenergy
      << "\t" << Benergy
      << "\t" << gas_energy << endl;
    my_file.close();
  }
}

void outputWrapperTXT::append_to_velocity_distribution(
  int cycle, int is, double maxVel, long long *VelocityDist)
{
  ofstream my_file(filename_ds.c_str(), fstream::app);
  my_file << cycle << "\t" << is << "\t" << maxVel;
  for (int i = 0; i < nbins; i++)
    my_file << "\t" << VelocityDist[i];
  my_file << endl;
  my_file.close();
}

// this assumes 4 species used in a certain way
//
void outputWrapperTXT::append_to_satellite_traces(Grid3DCU *grid,
    const_arr3_double Bx, const_arr3_double By, const_arr3_double Bz,
    const_arr3_double Ex, const_arr3_double Ey, const_arr3_double Ez,
    const_arr3_double Jxs,const_arr3_double Jys,const_arr3_double Jzs,
    const_arr3_double rhons)
{
  assert(do_write_virtual_satellite_traces());
  ofstream my_file(filename_cqsat.c_str(), fstream::app);
  const int nx0 = grid->get_nxc_r();
  const int ny0 = grid->get_nyc_r();
  const int nz0 = grid->get_nzc_r();
  const int nsat = get_number_of_satellites_per_direction();
  for (int isat = 0; isat < nsat; isat++)
  for (int jsat = 0; jsat < nsat; jsat++)
  for (int ksat = 0; ksat < nsat; ksat++)
  {
    const int i1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
    const int i2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
    const int i3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
    my_file << Bx.get(i1, i2, i3) << "\t"
            << By.get(i1, i2, i3) << "\t"
            << Bz.get(i1, i2, i3) << "\t";
    my_file << Ex.get(i1, i2, i3) << "\t"
            << Ey.get(i1, i2, i3) << "\t"
            << Ez.get(i1, i2, i3) << "\t";
    my_file << Jxs(0, i1, i2, i3) + Jxs(2, i1, i2, i3) << "\t"
            << Jys(0, i1, i2, i3) + Jys(2, i1, i2, i3) << "\t"
            << Jzs(0, i1, i2, i3) + Jzs(2, i1, i2, i3) << "\t";
    my_file << Jxs(1, i1, i2, i3) + Jxs(3, i1, i2, i3) << "\t"
            << Jys(1, i1, i2, i3) + Jys(3, i1, i2, i3) << "\t"
            << Jzs(1, i1, i2, i3) + Jzs(3, i1, i2, i3) << "\t";
    my_file << rhons(0, i1, i2, i3) + rhons(2, i1, i2, i3) << "\t";
    my_file << rhons(1, i1, i2, i3) + rhons(3, i1, i2, i3) << "\t";
  }
  my_file << endl;
  my_file.close();
}
