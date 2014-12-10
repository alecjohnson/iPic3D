#ifndef _IPIC3D_H_
#define _IPIC3D_H_

#include "ipic_fwd.h"
#include "assert.h"
#include "MIsolver.h"

namespace iPic3D
{
  // solver that chooses initial and boundary
  // conditions based on configuration
  // (move this into apps directory).
  //
  class c_Solver: public MIsolver
  {
   public:
    c_Solver(int argc, const char **argv) :
      MIsolver(argc, argv)
    {}
    // this class overrides this method
    virtual void set_initial_conditions();
  };
}

#endif
