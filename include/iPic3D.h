#ifndef _IPIC3D_H_
#define _IPIC3D_H_
// This file declares the following classes:
// 
//   MIsolver
//   c_Solver

#include "ipic_fwd.h"
#include "assert.h"

namespace iPic3D
{
  class c_Solver: public MIsolver
  {
    // this class overrides this method
    virtual void set_initial_conditions();
  };
}

#endif
