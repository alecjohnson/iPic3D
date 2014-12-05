
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

/* Interpolation smoothing: Smoothing (vector must already have ghost cells)
   value: if 1 nothing is done; else smoothing is done and value is ignored.
   type = 0 --> center based vector ;
   type = 1 --> node based vector ; */
void smooth(double value, arr3_double vector, int type,
  const VirtualTopology3D *vct,
  const Grid *grid)
{
  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {

    if (value != 1.0) {
      double alpha;
      int nx, ny, nz;
      switch (type) {
        case (0):
          nx = grid->getNXC();
          ny = grid->getNYC();
          nz = grid->getNZC();
          communicateCenterBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);

          break;
        case (1):
          nx = grid->getNXN();
          ny = grid->getNYN();
          nz = grid->getNZN();
          communicateNodeBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);
          break;
      }
      double ***temp = newArr3(double, nx, ny, nz);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      for (int i = 1; i < nx - 1; i++)
      for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++)
      {
          temp[i][j][k] = value * vector[i][j][k] + alpha *
             (vector[i - 1][j][k]
            + vector[i + 1][j][k]
            + vector[i][j - 1][k]
            + vector[i][j + 1][k]
            + vector[i][j][k - 1]
            + vector[i][j][k + 1]);
      }
      for (int i = 1; i < nx - 1; i++)
      for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++)
      {
          vector[i][j][k] = temp[i][j][k];
      }
      delArr3(temp, nx, ny);
    }
  }
}
/* Interpolation smoothing of 3-component node-based vector that has ghost nodes
   value: if 1 nothing is done; else smoothing is done and value is ignored.
   (could just call smooth three times, once for each component, passing in
   boundary conditions for electric field...)
 */
void smoothE(double value,
  arr3_double Ex;
  arr3_double Ey;
  arr3_double Ez;
  const Collective *col,
  const VirtualTopology3D *vct)
{
  const int nxn = grid->get_nxn();
  const int nyn = grid->get_nyn();
  const int nzn = grid->get_nzn();
  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Exth
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          temp[i][j][k] = value * Ex[i][j][k] + alpha *
              (Ex[i - 1][j][k]
             + Ex[i + 1][j][k]
             + Ex[i][j - 1][k]
             + Ex[i][j + 1][k]
             + Ex[i][j][k - 1]
             + Ex[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          Ex[i][j][k] = temp[i][j][k];
      // Eyth
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          temp[i][j][k] = value * Ey[i][j][k] + alpha *
              (Ey[i - 1][j][k]
             + Ey[i + 1][j][k]
             + Ey[i][j - 1][k]
             + Ey[i][j + 1][k]
             + Ey[i][j][k - 1]
             + Ey[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          Ey[i][j][k] = temp[i][j][k];
      // Ezth
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          temp[i][j][k] = value * Ez[i][j][k] + alpha *
              (Ez[i - 1][j][k]
             + Ez[i + 1][j][k]
             + Ez[i][j - 1][k]
             + Ez[i][j + 1][k]
             + Ez[i][j][k - 1]
             + Ez[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++)
          Ez[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
}

// === Section: MIsolver_routines ===

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
}

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

MIsolver::MIsolver(int argc, char **argv)
: setting(argc, argv),
  iMoments(0),
  EMf(0),
  kinetics(0),
  fieldForPcls(0),
  my_clock(0),
{
  ns = col->getNs(); // get the number of particle species involved in simulation

  #ifdef BATSRUS
  // set index offset for each processor
  setGlobalStartIndex(vct);
  #endif
}

int MIsolver::Init()
{
  #if defined(__MIC__)
  assert_eq(DVECWIDTH,8);
  #endif

  iMoments = new Imoments(setting);
  EMf = new EMfields3D(col, grid, vct);
  kinetics = new Kinetics(setting);

  set_initial_conditions();

  initialize_output();

  my_clock = new Timing(vct->get_rank());

  return 0;
}

Imoments::compute_from_primitive_moments(const Pmoments& pMoments,
  arr3_double Bx, arr3_double By, arr3_double Bz)
{
  //setZero();
  // sum all over the species
  //iMoments.sumOverSpecies();

  // copy the charge densities from the primitive moments
  for(int is=0;is<ns;is++)
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    rhons[is][i][j][k] = pMoments->get_rhons()[is][i][j][k];
  }

  // modify the charge density used to drive the electromagnetic
  // field based on non-kinetic influences.
  {
    // Fill with constant charge the planet
    if (setting.get_col().getCase()=="Dipole") {
      ConstantChargePlanet();
    }
    // Set a constant charge in the OpenBC boundaries
    ConstantChargeOpenBC();
  }

  calculateJhat(Bx, By, Bz);
}

void MIsolver::computeMoments()
{
  kinetics->compute_pMoments();

  // iMoments are the moments needed to drive the field
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
  iMoments.compute_from_primitive_moments(pMoments);

  // calculate densities on centers from nodes
  // (only needed for the Poisson
  // solve in the divergence cleaning, so probably
  // it would be better to do it there).
  EMf->interpDensitiesN2C();
  // calculate the hat quantities for the implicit method
  EMf->calculateRhoHat();
}
void Kinetics::compute_pMoments()
{
  timeTasks_set_main_task(TimeTasks::MOMENTS);

  // to be performed on booster
  {
    pad_particle_capacities();
    #pragma omp parallel
    for (int is = 0; is < ns; is++)
    {
      pMoments.accumulateMoments(part[is]);
      //#pragma omp master
      //[...send accumulated moments to cluster...]
    }
  }

  // to be performed on cluster
  for (int is = 0; is < ns; is++)
  {
    //[...wait for communicated moments to be received from booster...]
    pMoments.communicateGhostP2G(is);
  }
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
  EMf->set_fieldForPcls(kinetics.fetch_fieldForPcls());
}

//! MAXWELL SOLVER for Efield
void MIsolver::advanceEfield() {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  if(I_am_field_solver())
  {
    // define fields to use at open boundaries
    // based on magnetic field and Ohm's law
    EMf->updateInfoFields();
    // advance the E field
    EMf->calculateE();
  }
  send_field_to_kinetic_solver();
}

//! update Bfield (assuming Eth has already been calculated)
//  B^{n+1} = B^n - curl(Eth)
void MIsolver::advanceBfield() {
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
void MIsolver::moveParticles()
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
  writeRESTART(vct->get_rank(), cycle, ns, vct, col, grid, EMf, part, 0);
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
  // Velocity distribution
  if(!do_write_velocity_distribution()) return;
  //if(cycle % col->getVelocityDistributionOutputCycle() == 0)
  {
    for (int is = 0; is < ns; is++) {
      double maxVel = part[is].getMaxVelocity();
      const int nbins = OutputWrapperTXT::get_number_of_distribution_bins();
      long long *VelocityDist = part[is].getVelocityDistribution(nbins, maxVel);
      if (vct->is_rank0())
      {
        outputWrapperTXT->append_to_velocity_distribution(cycle, is, maxVel, VelocityDist);
      }
      delete [] VelocityDist;
    }
  }
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
    writeRESTART(vct->get_rank(), col->getNcycles() + col->getLast_cycle(), ns, vct, col, grid, EMf, part, 0);
    #endif
  }

  // stop profiling
  my_clock->stopTiming();
}

