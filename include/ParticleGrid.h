#ifndef _ParticleGrid_h_
#define _ParticleGrid_h_

#include "Grid.h"

class ParticleGrid : public Grid3DCU
{
 public:
  // inline methods to calculate mesh cell and weights.
  static void get_weights(double weights[8],
    double w0x, double w0y, double w0z,
    double w1x, double w1y, double w1z)
  {
    // which of the following is faster?
    //
    // this:
    //
    //const double weight00 = w0x*w0y;
    //const double weight01 = w0x*w1y;
    //const double weight10 = w1x*w0y;
    //const double weight11 = w1x*w1y;
    //weights[0] = weight00*w0z; // weight000
    //weights[1] = weight00*w1z; // weight001
    //weights[2] = weight01*w0z; // weight010
    //weights[3] = weight01*w1z; // weight011
    //weights[4] = weight10*w0z; // weight100
    //weights[5] = weight10*w1z; // weight101
    //weights[6] = weight11*w0z; // weight110
    //weights[7] = weight11*w1z; // weight111
    //
    // or this:
    //
    weights[0] = w0x*w0y*w0z; // weight000
    weights[1] = w0x*w0y*w1z; // weight001
    weights[2] = w0x*w1y*w0z; // weight010
    weights[3] = w0x*w1y*w1z; // weight011
    weights[4] = w1x*w0y*w0z; // weight100
    weights[5] = w1x*w0y*w1z; // weight101
    weights[6] = w1x*w1y*w0z; // weight110
    weights[7] = w1x*w1y*w1z; // weight111
  }
  void get_cell_coordinates(
    int& cx, int& cy, int& cz,
    double xpos, double ypos, double zpos)
  {
      // xStart marks start of domain excluding ghosts
      const double rel_xpos = xpos - xStart;
      const double rel_ypos = ypos - yStart;
      const double rel_zpos = zpos - zStart;
      // cell position minus 1 (due to ghost cells)
      const double cxm1_pos = rel_xpos * invdx;
      const double cym1_pos = rel_ypos * invdy;
      const double czm1_pos = rel_zpos * invdz;
      cx = 1 + int(floor(cxm1_pos));
      cy = 1 + int(floor(cym1_pos));
      cz = 1 + int(floor(czm1_pos));
  }
  void make_cell_coordinates_safe(int& cx, int& cy, int& cz)
  {
    // if the cell is outside the domain, then treat it as
    // in the nearest ghost cell.
    //
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    if (cz < 0) cz = 0;
    if (cx > cxlast) cx = cxlast; //nxc-1;
    if (cy > cylast) cy = cylast; //nyc-1;
    if (cz > czlast) cz = czlast; //nzc-1;
  }
  void get_safe_cell_coordinates(
    int& cx, int& cy, int& cz,
    double x, double y, double z)
  {
    get_cell_coordinates(cx,cy,cz,x,y,z);
    make_cell_coordinates_safe(cx,cy,cz);
  }
  void get_safe_cell_and_weights(
    double xpos, double ypos, double zpos,
    int &cx, int& cy, int& cz,
    double weights[8])
  {
    //convert_xpos_to_cxpos(xpos,ypos,zpos,cxpos,cypos,czpos);
    // xStart marks start of domain excluding ghosts
    const double rel_xpos = xpos - xStart;
    const double rel_ypos = ypos - yStart;
    const double rel_zpos = zpos - zStart;
    // cell position minus 1 (due to ghost cells)
    const double cxm1_pos = rel_xpos * invdx;
    const double cym1_pos = rel_ypos * invdy;
    const double czm1_pos = rel_zpos * invdz;
    //
    cx = 1 + int(floor(cxm1_pos));
    cy = 1 + int(floor(cym1_pos));
    cz = 1 + int(floor(czm1_pos));
  
    make_cell_coordinates_safe(cx,cy,cz);
  
    // fraction of the distance from the right of the cell
    const double w1x = cx - cxm1_pos;
    const double w1y = cy - cym1_pos;
    const double w1z = cz - czm1_pos;
    // fraction of distance from the left
    const double w0x = 1.-w1x;
    const double w0y = 1.-w1y;
    const double w0z = 1.-w1z;

    get_weights(weights, w0x, w0y, w0z, w1x, w1y, w1z);
  }
  void get_safe_cell_and_weights(double xpos[3], int cx[3], double weights[8])
  {
    get_safe_cell_and_weights(xpos[0],xpos[1],xpos[2],cx[0],cx[1],cx[2],weights);
  }

  // methods used in sorting particles
  //
  //int get_num_cells_per_thread() { return num_cells_per_thread; }
  // would it be faster to look this up in an array?
  int get_thread_for_cell(int idx) {
    return idx/num_cells_per_thread;
  }
  //
  // calculate cell of particle
  //
  // (These need to be inline for fast processing.)
  // get one-dimensional index of cell at position
  int get_cell_idx_for_coordinates(int cx, int cy, int cz)
  {
    return cz+nzc*(cy+nyc*cx);
  }
  int get_safe_cell_idx(double x, double y, double z)
  {
    int cx, cy, cz;
    get_safe_cell_coordinates(cx,cy,cz,x,y,z);
    return get_cell_idx_for_coordinates(cx, cy, cz);
  }
  /*! version that assumes particle is in domain */
  void get_cell_for_pos_in_domain(
    int& cx, int& cy, int& cz, 
    pfloat xpos, pfloat ypos, pfloat zpos)
  {
    // xstart is left edge of domain excluding ghost cells
    // cx=0 for ghost cell layer.
    cx = 1 + int(floor((xpos - xStart) * invdx));
    cy = 1 + int(floor((ypos - yStart) * invdy));
    cz = 1 + int(floor((zpos - zStart) * invdz));
    //
    assert_le(0,cx);
    assert_le(0,cy);
    assert_le(0,cz);
    assert_le(cx,nxc);
    assert_le(cy,nyc);
    assert_le(cz,nzc);
  }
};

#endif
