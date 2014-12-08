
// === Section: Kinetics_routines ===

Kinetics::~Kinetics()
{
  // delete particles
  //
  if(part)
  {
    for (int i = 0; i < ns; i++)
    {
      // placement delete
      part[i].~Particles3D();
    }
    free(part);
  }
}

Kinetics::Kinetics(const Setting& setting_)
: setting(setting_)
{
  // Allocation of particles
  // part = new Particles3D[ns];
  part = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    // placement new
    new(&part[i]) Particles3D(i,setting);
  }
}

/*  -------------- */
/*!  Particle mover */
/*  -------------- */
bool Kinetics::moveParticles(const_arr4_double fieldForPcls)
{
  // move all species of particles
  {
    timeTasks_set_main_task(TimeTasks::PARTICLES);

    pad_particle_capacities();

    #pragma omp parallel
    {
    for (int i = 0; i < ns; i++)  // move each species
    {
      // #pragma omp task inout(part[i]) in(grid) target_device(booster)
      //
      // use the Predictor Corrector scheme to move particles
      switch(Parameters::get_MOVER_TYPE())
      {
        case Parameters::SoA:
          part[i].mover_PC(fieldForPcls);
          break;
        //case Parameters::SoA_vec_resort:
        //  part[i].mover_PC_vectorized(fieldForPcls);
        //  break;
        case Parameters::AoS:
          part[i].mover_PC_AoS(fieldForPcls);
          break;
        case Parameters::AoSintr:
          part[i].mover_PC_AoS_vec_intr(fieldForPcls);
          break;
        case Parameters::AoSvec:
          part[i].mover_PC_AoS_vec(fieldForPcls);
          break;
        //case Parameters::AoS_vec_onesort:
        //  part[i].mover_PC_AoS_vec_onesort(fieldForPcls);
        //  break;
        default:
          unsupported_value_error(Parameters::get_MOVER_TYPE());
      }
      // overlap initial communication of electrons with moving of ions
      #pragma omp master
      part[i].separate_and_send_particles();
    }
    }
    for (int i = 0; i < ns; i++)  // communicate each species
    {
      //part[i].communicate_particles();
      part[i].recommunicate_particles_until_done(1);
    }
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */

  for (int i=0; i < ns; i++) {
    if (col->getRHOinject(i)>0.0)
      part[i].repopulate_particles();
  }

  /* --------------------------------------- */
  /* Remove particles from depopulation area */
  /* --------------------------------------- */

  if (col->getCase()=="Dipole") {
    for (int i=0; i < ns; i++)
    {
      double Qremoved = part[i].deleteParticlesInsideSphere( col->getL_square(),
        col->getx_center(),col->gety_center(),col->getz_center());
    }
  }
  return (false);
}

void Kinetics::sortParticles() {
  eprintf("check how this is being used.");
  for(int species_idx=0; species_idx<ns; species_idx++)
    part[species_idx].sort_particles();
}

// convert particle to struct of arrays (assumed by I/O)
void Kinetics::convertParticlesToSoA()
{
  eprintf("check how this is being used.");
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSoA();
}

// convert particle to array of structs (used in computing)
void Kinetics::convertParticlesToAoS()
{
  eprintf("check how this is being used.");
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToAoS();
}

// convert particle to array of structs (used in computing)
void Kinetics::convertParticlesToSynched()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSynched();
}

void Kinetics::compute_pMoments(Pmoments& pMoments)
{
  timeTasks_set_main_task(TimeTasks::MOMENTS);
  // to be performed on booster
  fetch_Pmoments().accumulateMoments(part);
  // to be performed on cluster
  fetch_Pmoments().communicateGhostP2G();
}
