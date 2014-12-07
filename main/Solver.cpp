
#include "mpi.h"
#ifndef NO_HDF5
#include "Restart3D.h"
#endif
#include "MPIdata.h"
#include "iPic3D.h"
#include "TimeTasks.h"
#include "ipicdefs.h"
#include "debug.h"
#include "Parameters.h"
#include "ompdefs.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "Moments.h"
#include "Timing.h"
//
#ifndef NO_HDF5
#include "ParallelIO.h"
#include "WriteOutputParallel.h"
#include "OutputWrapperFPP.h"
#endif
//
#include <iostream>
#include <sstream>

using namespace iPic3D;

// === Section: MIsolver_routines ===

// todo: separate out this problem-specific code
//
void c_Solver::set_initial_conditions()
{
  if (col->getCase()=="Dipole")         initDipole();
  // cases prior to this point are responsible for their own restarts
  else if(col->getRestart_status())     init_from_restart();
  // cases that use rhon to set rhoc
  else if (col->getCase()=="restart")   init_from_restart();
  else if (col->getCase()=="GEM")       initGEM();
  else if (col->getCase()=="GEMnoPert") initGEMnoPert();
  else if (col->getCase()=="ForceFree") initForceFree();
  else if (col->getCase()=="GEM")       initGEM();
  // cases that use rhoc to set rhon
#ifdef BATSRUS
  else if (col->getCase()=="BATSRUS")   initBATSRUS();
#endif
  else if (col->getCase()=="RandomCase")initRandomField();
  else {
    eprintf("Case=%s in inputfile is not supported.", col->getCase());
  }
  //if (myrank==0) {
  //  cout << "Case is " << col->getCase() <<"\n";
  //  cout <<"total # of particle per cell is " << col->getNpcel(0) << "\n";
  //}

  // OpenBC
  EMf->updateInfoFields();

  // Allocation of particles
  // part = new Particles3D[ns];
  part = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    new(&part[i]) Particles3D(i,col,vct,grid);
  }
  //particleSolver = new ParticleSolver(col,vct,grid);

  // Initial Condition for PARTICLES if you are not starting from RESTART
  if(col->getRestart_status() == 0)
  {
    // wave = new Planewave(col, EMf, grid, vct);
    // wave->Wave_Rotated(part); // Single Plane Wave
    for (int i = 0; i < ns; i++)
    {
      if      (col->getCase()=="ForceFree") part[i].force_free(EMf);
#ifdef BATSRUS
      else if (col->getCase()=="BATSRUS")   part[i].MaxwellianFromFluid(col,i);
#endif
      else                                  part[i].maxwellian(EMf);
      part[i].reserve_remaining_particle_IDs();
    }
  }
}

void MIsolver::initialize_output()
{
  if (col->getWriteMethod() == "shdf5")
  {
    #ifndef NO_HDF5
    outputWrapperFPP = new OutputWrapperFPP(col,vct,grid);
    fetch_outputWrapperFPP().init_output_files(EMf,part);
    #endif
  }
  outputWrapperTXT = new OutputWrapperTXT(col,vct);
}

MIsolver::~MIsolver()
{
  delete pMoments;
  delete EMf; // field
  delete kinetics;
  delete fieldForPcls;
  #ifndef NO_HDF5
  delete outputWrapperFPP;
  delete outputWrapperTXT;
  #endif
  delete my_clock;
  delete &setting;
}


// MIsolver::MIsolver(int argc, char **argv)
// : setting(argc, argv),
//   iMoments(*new Imoments(setting)),
//   EMf(*new EMfields3D(setting, iMoments)),
//   pMoments(*new Pmoments(setting)),
//   kinetics(*new Kinetics(setting)),
//   fieldForPcls(*new array4_double(
//     setting.grid().get_nxn(),
//     setting.grid().get_nyn(),
//     setting.grid().get_nzn(),
//     2*DFIELD_3or4)),
//   my_clock(0)
MIsolver::MIsolver(int argc, char **argv)
: setting(*new Setting(argc, argv)),
  iMoments(0),
  EMf(0),
  pMoments(0),
  kinetics(0),
  fieldForPcls(0),
  my_clock(0),
{
  #ifdef BATSRUS
  // set index offset for each processor
  setGlobalStartIndex(vct);
  #endif

  #if defined(__MIC__)
  assert_eq(DVECWIDTH,8);
  #endif

  iMoments = new Imoments(setting);
  EMf = new EMfields3D(setting, iMoments);
  kinetics = new Kinetics(setting);
  pMoments = new Pmoments(setting);

  set_initial_conditions();

  initialize_output();

  my_clock = new Timing(vct->get_rank());

  return 0;
}

void MIsolver::compute_moments()
{
  get_kinetics().compute_pMoments(pMoments);

  // iMoments are the moments needed to drive the implicit field
  // solver, which are separated out from EMfields3D in case
  // we someday wish to communicate these 3+num_species
  // moments instead of the 10*num_species primitive
  // moments.  Note that computing implicit moments would
  // require communicating boundary data for 6*num_species
  // pressure tensor components.  Since the current
  // implmentation of boundary communication has many
  // synchronizations, for the present we prefer to avoid
  // all exchange of boundary data on the booster.
  //
  // by default we smooth after applying magnetic field
  const_arr3_double Bx = EMf->get_Bx_smooth();
  const_arr3_double By = EMf->get_By_smooth();
  const_arr3_double Bz = EMf->get_Bz_smooth();
  // smooth before applying magnetic field
  if(Parameters::use_perfect_smoothing())
  {
    Bx = EMf->get_Bx_tot();
    By = EMf->get_By_tot();
    Bz = EMf->get_Bz_tot();
  }
  iMoments->compute_from_primitive_moments(pMoments, Bx, By, Bz);

  // could wait to do this until field solve
  //
  // calculate the hat quantities for the implicit method
  EMf->calculateRhoHat();
}

// This method should send or receive field
// depending on whether the process is
// a cluster process or a booster process.
// The proper way to do this (rather than
// using if statements) would be to make
// this a virtual method in an MIsolver_base
// class from which MIsolver, MIfieldSolver,
// and MIkineticSolver all inherit.
//
// Note that if will be more efficient to send
// the SoA fields and then transpose them than
// to communicate the AoS fields, because
// data is padded by 33% in the SoA fields
// representation for alignment purposes,
// Moreover, the magnetic field could be sent
// much earlier than the electric field
// to hide bandwidth limitations (in which
// case this method might more accurately be named
// "finish_transferring_field_to_kinetic_solver()".
// The expense of transposing the field data
// on the Booster should be minor compared to 
// the reduction in communication expense.
//
void MIsolver::send_field_to_kinetic_solver()
{
  set_fieldForPcls();
}
void MIsolver::set_fieldForPcls()
{
  //EMf->set_fieldForPcls(fetch_fieldForPcls());
  #pragma omp parallel for collapse(2)
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    fieldForPcls[i][j][k][0] = Bx_smooth[i][j][k];
    fieldForPcls[i][j][k][1] = By_smooth[i][j][k];
    fieldForPcls[i][j][k][2] = Bz_smooth[i][j][k];
    fieldForPcls[i][j][k][0+DFIELD_3or4] = Ex_smooth[i][j][k];
    fieldForPcls[i][j][k][1+DFIELD_3or4] = Ey_smooth[i][j][k];
    fieldForPcls[i][j][k][2+DFIELD_3or4] = Ez_smooth[i][j][k];
  }
}

//! MAXWELL SOLVER for Efield
void MIsolver::advance_Efield() 
{
  timeTasks_set_main_task(TimeTasks::FIELDS);
  if(I_am_field_solver())
  {
    // advance the E field
    EMf->calculateE(get_iMoments());
  }
  send_field_to_kinetic_solver();
}

//! update Bfield (assuming Eth has already been calculated)
//  B^{n+1} = B^n - curl(Eth)
void MIsolver::advance_Bfield() {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  timeTasks_set_task(TimeTasks::BFIELD); // subtask
  // calculate the B field
  EMf->advanceB();
  // coupling: send_Bfield_to_kinetics();
  // (This communication should be easily hidden,
  // since the solver does not need to receive Bfield
  // until just prior to the particle move.)
}

// this method should be a no-op on the cluster
void MIsolver::move_particles()
{
  //[...receive field from fieldsolver...]
  kinetics->moveParticles();
}

// --- section: output methods ---

void MIsolver::WriteRestart(int cycle)
{
  bool do_WriteRestart = (cycle % col->getRestartOutputCycle() == 0 && cycle != col->get_first_cycle());
  if(!do_WriteRestart)
    return;

  #ifndef NO_HDF5
  kinetics->convertParticlesToSynched(); // hack
  // write the RESTART file
  // without 0 add to restart file
  writeRESTART(vct->get_rank(), cycle, vct, col, grid, EMf, part, 0);
  #endif
}

// write the conserved quantities
void MIsolver::WriteConserved(int cycle) {
  if(cycle % col->getDiagnosticsOutputCycle() == 0)
  {
    double Eenergy = EMf->getEenergy();
    double Benergy = EMf->getBenergy();
    double gas_energy = 0.0;
    double bogus_momentum = 0.0;
    for (int is = 0; is < ns; is++) {
      gas_energy += part[is].getKe();
      bogus_momentum += part[is].getP();
    }
    double total_energy = gas_energy + Eenergy + Benergy;
    outputWrapperTXT->append_conserved_quantities(
      total_energy, bogus_momentum, Eenergy, Benergy, gas_energy);
  }
}

void MIsolver::WriteVelocityDistribution(int cycle)
{
  outputWrapperTXT->write_velocity_distribution(cycle, part);
}

// This seems to record values at a grid of sample points
//
void MIsolver::WriteVirtualSatelliteTraces()
{
  if(!do_write_virtual_satellite_traces()) return;

  assert_eq(ns,4);
  outputWrapperTXT->append_to_satellite_traces(grid,
    EMf->getBx(), EMf->getBy(), EMf->getBz(),
    EMf->getEx(), EMf->getEy(), EMf->getEz(),
    EMf->getJxs(), EMf->getJys(), EMf->getJzs(),
    EMf->get_rhons());
}

void MIsolver::WriteFields(int cycle) {
  if(col->field_output_is_off())
    return;
  #ifdef NO_HDF5
    eprintf("must compile with HDF5");
  #else
  if(cycle % (col->getFieldOutputCycle()) == 0 || cycle == col->get_first_cycle())
  {
    timeTasks_set_task(TimeTasks::WRITE_FIELDS);
    if (col->getWriteMethod() == "Parallel") {
        WriteOutputParallel(grid, EMf, col, vct, cycle);
    }
    else // OUTPUT to large file, called proc**
    {
        // Pressure tensor is available
        fetch_outputWrapperFPP().append_output(
          "Eall + Ball + rhos + Jsall + pressure", cycle);
    }
  }
  #endif
}

void MIsolver::WriteParticles(int cycle)
{
  if(col->particle_output_is_off())
    return;
  #ifdef NO_HDF5
    eprintf("NO_HDF5 requires OutputMethod=none")
  #else

  timeTasks_set_task(TimeTasks::WRITE_PARTICLES);

  // this is a hack
  kinetics->convertParticlesToSynched();

  if (col->getWriteMethod() == "Parallel")
  {
    dprintf("pretending to write particles (not yet implemented)");
  }
  else
  {
    fetch_outputWrapperFPP().append_output(
      "position + velocity + q ", cycle, 1);
  }
  #endif // NO_HDF5
}

// This needs to be separated into methods that save particles
// and methods that save field data
//
void MIsolver::WriteOutput(int cycle) {

  // The quickest things should be written first.

  WriteConserved(cycle);
  WriteVelocityDistribution(cycle);

  // mechanism to suppress output
  if(!Parameters::get_doWriteOutput())
    return;

  #ifndef NO_HDF5 // array output is only supported for HDF5
  // once phdf5 is properly supported,
  // let's change "Parallel" to mean "phdf5".
  if (col->getWriteMethod() == "H5hut"
   || col->getWriteMethod() == "Parallel")
  {
    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */

    if (!col->field_output_is_off() &&
        cycle%(col->getFieldOutputCycle())==0)
      WriteFieldsH5hut(ns, grid, EMf, col, vct, cycle);
    if (!col->particle_output_is_off() &&
        cycle%(col->getParticlesOutputCycle())==0)
      WritePartclH5hut(ns, grid, part, col, vct, cycle);
  }
  else if (col->getWriteMethod() == "phdf5")
  {
    /* -------------------------------------------- */
    /* Parallel output using basic hdf5 */
    /* -------------------------------------------- */

    if (!col->field_output_is_off() &&
      cycle%(col->getFieldOutputCycle())==0)
      WriteOutputParallel(grid, EMf, part, col, vct, cycle);
    if (!col->particle_output_is_off() &&
        cycle%(col->getParticlesOutputCycle())==0)
    {
      if(!MPIdata::get_rank())
        warning_printf("WriteParticlesParallel() is not yet implemented.");
      //WritePartclH5hut(ns, grid, part, col, vct, cycle);
    }
  }
  else if (col->getWriteMethod() == "shdf5")
  {
    // write fields-related data
    WriteFields(cycle);
    // This should be invoked by user if desired
    // by means of a callback mechanism.
    WriteVirtualSatelliteTraces();

    // write particles-related data
    //
    // this also writes field data...
    WriteRestart(cycle);
    WriteParticles(cycle);
  }
  else if (col->getWriteMethod() == "default")
  {
    if(col->getParticlesOutputCycle()==1)
    {
      warning_printf(
        "ParticlesOutputCycle=1 now means output particles with evey cycle.\n"
        "\tParticlesOutputCycle = 0 turns off particle output.");
    }
    eprintf("The new name for serial hdf5 output is shdf5.\n"
      "\tselect WriteMethod=shdf5.");
  }
  else
  {
    invalid_value_error(col->getWriteMethod().c_str());
  }
  #endif
}

void MIsolver::Finalize() {
  if (col->getCallFinalize())
  {
    #ifndef NO_HDF5
    get_kinetics().convertParticlesToSynched();
    writeRESTART(vct->get_rank(), col->getNcycles() + col->getLast_cycle(), vct, col, grid, EMf, part, 0);
    #endif
  }

  // stop profiling
  my_clock->stopTiming();
}

// Flow structure of application.
// C = "lives on Cluster"
// B = "lives on Booster"
//
// C    En          = advance(dt, En, Bn, Btot, Jhat, rhon)
// C    Bc          = advance(dt, Bc, En)
// C    Bn          = interpolate_to_nodes(Bc)
// C    Btot        = Bn + Bext
// C B  B_smooth    = smooth(Btot)
// C B  E_smooth    = smooth(En)
//   B  BEaos       = soa2aos(B_smooth,E_smooth)
//   B  particles   = advance(dt, particles, BEaos)
// C B  pMoments    = sumMoments(particles)
// C    Jhat_coarse = compute_Jhat(pMoments, B_smooth)
// C    Jhat        = smooth(compute_Jhat(pMoments, B_smooth)
// C    rhon        = smooth(pMoments.rhon)
//
// remarks:
//  
// * we could transfer Jhat_coarse rather than pMoments
// * a more correct way to compute Jhat is:
//   C    Jhat        = compute_Jhat(smooth(Pmoments), Btot)
// * Btot is used to compute MUdot in En.advance().
// * Btot or B_smooth is used in Jhat.smooth() to compute PIdot.
// * Jhat could be transferred from B to C rather than Pmoments.
// * smoothing requires exchange of boundary data.
// * computing Jhat requires exchange of boundary data.
//
// Class organization:
// MIsolver:
//   EMfields3D:
//     En, Bc, Bn, Btot
//   B_smooth
//   E_smooth
//   BEaos
//   pMoments
//   particles
// FieldSolver:
//   EMfields3D:
//     En, Bc, Bn, Btot
//   B_smooth
//   E_smooth
//   pMoments
//   Jhat_coarse
//   Jhat
//   rhon
// KineticSolver:
//   B_smooth
//   E_smooth
//   BEaos
//   particles
//   pMoments
//
void MIsolver::run()
{
  timeTasks.resetCycle();
  compute_moments();
  for (int i = FirstCycle(); i <= FinalCycle(); i++)
  {
    if (is_rank0())
      printf(" ======= Cycle %d ======= \n",i);

    timeTasks.resetCycle();
    advance_Efield();
    move_particles();
    advance_Bfield();
    compute_moments();
    WriteOutput(i);
    // print out total time for all tasks
    timeTasks.print_cycle_times(i);
  }

  Finalize();
}
