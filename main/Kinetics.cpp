#include "Kinetics.h"
#include "SpeciesMoms.h"
#include "Particles3D.h"
#include "Setting.h"
#include "Collective.h"
#include "TimeTasks.h"
#include <new> // needed for placement new

Kinetics::~Kinetics()
{
  // delete particles
  //
  if(speciesPcls)
  {
    for (int is = 0; is < ns; is++)
    {
      // placement delete
      speciesPcls[is].~Particles3D();
    }
    free(speciesPcls);
  }
}

Kinetics::Kinetics(const Setting& setting_)
: setting(setting_), ns(setting.col().getNs())
{
  // Allocation of particles
  // speciesPcls = new Particles3D[ns];
  speciesPcls = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int is = 0; is < ns; is++)
  {
    // placement new
    new(&speciesPcls[is]) Particles3D(is,setting);
  }
}

Particles3D& Kinetics::fetch_pcls(int is)
{
  return speciesPcls[is];
}

/*  -------------- */
/*!  Particle mover */
/*  -------------- */
bool Kinetics::moveParticles(const_arr4_double fieldForPcls)
{
  const Collective& col = setting.col();
  // move all species of particles
  {
    timeTasks_set_main_task(TimeTasks::PARTICLES);
    for (int is = 0; is < ns; is++)
      speciesPcls[is].pad_capacities();
    #pragma omp parallel
    {
    for (int is = 0; is < ns; is++)  // move each species
    {
      // #pragma omp task inout(speciesPcls[is]) in(grid) target_device(booster)
      //
      // use the Predictor Corrector scheme to move particles
      switch(Parameters::get_MOVER_TYPE())
      {
        case Parameters::SoA:
          speciesPcls[is].mover_PC(fieldForPcls);
          break;
        //case Parameters::SoA_vec_resort:
        //  speciesPcls[is].mover_PC_vectorized(fieldForPcls);
        //  break;
        case Parameters::AoS:
          speciesPcls[is].mover_PC_AoS(fieldForPcls);
          break;
        case Parameters::AoSintr:
          speciesPcls[is].mover_PC_AoS_vec_intr(fieldForPcls);
          break;
        case Parameters::AoSvec:
          speciesPcls[is].mover_PC_AoS_vec(fieldForPcls);
          break;
        //case Parameters::AoS_vec_onesort:
        //  speciesPcls[is].mover_PC_AoS_vec_onesort(fieldForPcls);
        //  break;
        default:
          unsupported_value_error(Parameters::get_MOVER_TYPE());
      }
      // overlap initial communication of electrons with moving of ions
      #pragma omp master
      speciesPcls[is].separate_and_send_particles();
    }
    }
    for (int is = 0; is < ns; is++)  // communicate each species
    {
      //speciesPcls[is].communicate_particles();
      speciesPcls[is].recommunicate_particles_until_done(1);
    }
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */

  for (int is=0; is < ns; is++) {
    if (col.getRHOinject(is)>0.0)
      speciesPcls[is].repopulate_particles();
  }

  /* --------------------------------------- */
  /* Remove particles from depopulation area */
  /* --------------------------------------- */

  if (col.getCase()=="Dipole") {
    for (int is=0; is < ns; is++)
    {
      double Qremoved = speciesPcls[is].deleteParticlesInsideSphere(
        col.getL_square(),
        col.getx_center(),col.gety_center(),col.getz_center());
    }
  }
  return (false);
}

void Kinetics::sortParticles() {
  eprintf("check how this is being used.");
  for(int species_idx=0; species_idx<ns; species_idx++)
    speciesPcls[species_idx].sort_particles();
}

// convert particle to struct of arrays (assumed by I/O)
void Kinetics::convertParticlesToSoA()
{
  eprintf("check how this is being used.");
  for (int is = 0; is < ns; is++)
    speciesPcls[is].convertParticlesToSoA();
}

// convert particle to array of structs (used in computing)
void Kinetics::convertParticlesToAoS()
{
  eprintf("check how this is being used.");
  for (int is = 0; is < ns; is++)
    speciesPcls[is].convertParticlesToAoS();
}

// convert particle to array of structs (used in computing)
void Kinetics::convertParticlesToSynched()
{
  for (int is = 0; is < ns; is++)
    speciesPcls[is].convertParticlesToSynched();
}

void Kinetics::reserve_remaining_particle_IDs()
{
  for (int is = 0; is < ns; is++)
    speciesPcls[is].reserve_remaining_particle_IDs();
}

