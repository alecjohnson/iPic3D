#ifndef Setting_h
#define Setting_h

class Collective;
class VCtopology3D;
class Grid3DCU;

// This class wraps up the background environment of the solver,
// including configuration parameters, other processes, and the grid.
//
class Setting
{
  Collective    *_col;  // parameter access
  VCtopology3D  *_vct;  // MPI process cartesian topology
  Grid3DCU      *_grid; // mesh geometry

 public:
  Setting(int argc, const char** argv);
  ~Setting();

  //const Collective&   get_col ()const{return *col; }
  //const VCtopology3D& get_vct ()const{return *vct; }
  //const Grid3DCU&     get_grid()const{return *grid;}
  const Collective&   col ()const{return *_col; }
  const VCtopology3D& vct ()const{return *_vct; }
  const Grid3DCU&     grid()const{return *_grid;}
};

#endif // Setting_h
