#ifndef MImoments_H
#define MImoments_H
#include "ipicarrays.h"
#include "Setting.h"

class Particles3Dcomm;
class SpeciesMoms;

// moments that must be communicated from the kinetic solver
// to the field solver of the moment-implicit method.
class MImoments
{
  private: // parameters
    // access to global data
    const Setting& setting;
    // dimensions of arrays
    const int nxn, nyn, nzn;
    // number of species
    const int ns;
  private: // moments
    // densities of each species
    vector_array3_double rhons;

    // implicit current density defined on nodes
    //
    array3_double Jxh;
    array3_double Jyh;
    array3_double Jzh;

  public: // accessors
    const vector_array3_double& get_rhons()const{return rhons;}
    const_arr3_double get_Jxh()const{return Jxh;}
    const_arr3_double get_Jyh()const{return Jyh;}
    const_arr3_double get_Jzh()const{return Jzh;}

  public:
    MImoments(const Setting& setting_);
    void compute_from_speciesMoms(const SpeciesMoms& speciesMoms,
      const_arr3_double Bx, const_arr3_double By, const_arr3_double Bz);
  private:
    void ConstantChargeOpenBCv2();
    void ConstantChargeOpenBC();
    void ConstantChargePlanet();
};

#endif
