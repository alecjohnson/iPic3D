#ifndef PclBins_h
#define PclBins_h
#include "aligned_vector.h"
#include "Alloc.h"

class Grid3DCU;

class PclBins
{
  int num_bins[3] ALLOC_ALIGNED;
  int last_idx[3] ALLOC_ALIGNED;
  int cells_per_bin[3] ALLOC_ALIGNED;
  int beg_pos[3] ALLOC_ALIGNED;
  int invdpos[3] ALLOC_ALIGNED;
  int nbins;
  array3<vector_SpeciesParticle> pclBins;
 public:
  PclBins(
    Grid3DCU* grid_,
    int num_bins_[3])
  
 public: // accessors
  int get_nbins(){return nbins;}
  vector_SpeciesParticle& fetch(int xi, int yi, int zi)
  {
    return pclBins.fetch(xi,yi,zi);
  }
  vector_SpeciesParticle& fetch(int bi)
  {
    return pclBins.vfetch(bi);
  }
 public: // calculated access
  vector_SpeciesParticle& fetch(double x, double y, double z)
  {
    // will the alignment help vectorization of this
    // performance-sensitive calculation (for sorting
    // particles)?
    int bin[3] ALLOC_ALIGNED;
    for(int i=0;i<3;i++)
    {
      bin[i] = floor((pos[i]-beg[i])*invdpos[i]);
      // clamp bounds
      bin[i] = std::max(bin[i],0);
      bin[i] = std::min(bin[i],last_idx[i]);
    }
    return pclBins.fetch(bin[0],bin[1],bin[2]);
  }
  // return total number of particles in all bins
  int ret_nop()
  {
    int nop=0;
    for(int bi=0;bi<nbins;bi++)
      nop += pclBins.vfetch(bi).size();
    return nop;
  }
};

// inline methods that should be moved into PclBins.cpp

inline void PclBins::clear()
{
  vector_SpeciesParticle* pclBins_arr = pclBins.fetch_arr();
  for(int i=0;i<nbins;i++)
    pclBins_arr[i].clear();
}

inline PclBins::PclBins(
  Grid3DCU* grid_,
  int num_bins_[3])
:
  pclBins(num_bins[0], num_bins[1], num_bins[2]);
{
  num_bins[0] = num_bins_[0];
  num_bins[1] = num_bins_[1];
  num_bins[2] = num_bins_[2];
  nbins= num_bins[0]*num_bins[1]*num_bins[2];
  beg_pos[0] = grid_->get_beg_pos()[0];
  beg_pos[1] = grid_->get_beg_pos()[1];
  beg_pos[2] = grid_->get_beg_pos()[2];
  cells_per_bin[0] = ceiling_of_ratio(grid_->get_nxc_r(),num_bins[0]);
  cells_per_bin[1] = ceiling_of_ratio(grid_->get_nyc_r(),num_bins[1]);
  cells_per_bin[2] = ceiling_of_ratio(grid_->get_nzc_r(),num_bins[2]);
  invdpos[0] = grid_->get_invdx()*cells_per_bin[0];
  invdpos[1] = grid_->get_invdy()*cells_per_bin[1];
  invdpos[2] = grid_->get_invdz()*cells_per_bin[2];
  last_idx[0] = num_bin[0]-1;
  last_idx[1] = num_bin[1]-1;
  last_idx[2] = num_bin[2]-1;
}

#endif // PclBins_h
