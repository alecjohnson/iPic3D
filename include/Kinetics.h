/***************************
  Kinetics: complements the FieldSolver.
  Now contains list of ParticleSolver instances,
  one for each species.  I plan to extend this
  to include fluid species as well.  
  In the future I want to support Vlasov species
  and combination of kinetic models with fluid models
  that use kinetic closure.
 *************************** */

#ifndef _IPIC3D_H_
#define _IPIC3D_H_

class Timing;

#include "ipic_fwd.h"
#include "assert.h"

// I called this "Kinetics" rather than "Particles"
// because I hope to support evolving a fluid with
// kinetic closure.
class Kinetics
{
  private:
    // References to outside data
    const Setting&setting;
    // variables that this class is responsible for
    Particles3D  *part;
    //
    // put these in MIsolver
    //
    //SpeciesMoms      speciesMoms;
    //array4_double fieldForPcls; // rename BEaos?
    // Electric field component used to move particles
    // organized in AoS format for rapid random access in particle mover.
    //array4_double fieldForPcls;

  private: // accessors

  public:
    ~Kinetics();
    Kinetics(const Setting& setting_);
    // get hatted moments for field solver
    void compute_speciesMoms(SpeciesMoms& speciesMoms)const;
    //void calculateMoments();
    // advance particles in response to electromagnetic field
    bool moveParticles(const_arr4_double fieldForPcls);
    // update magnetic field using Eth as in CalculateB
    // B^{n+1} = B^n - curl(Eth)
    //void UpdateB();
    // if we use Eth to push the particles then we don't need this at all.
    // E^{n+1} = (1/theta)*Eth + (1/theta-1)*E^{n}
    //void UpdateE();
    //
  private:
    void pad_particle_capacities();
    void convertParticlesToSoA();
    void convertParticlesToAoS();
    void convertParticlesToSynched();
    void sortParticles();
};

#endif
