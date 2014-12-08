
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

MIsolver::~MIsolver()
{
  delete speciesMoms;
  delete miMoments;
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
//   miMoments(*new MImoments(setting)),
//   EMf(*new EMfields3D(setting, miMoments)),
//   speciesMoms(*new SpeciesMoms(setting)),
//   kinetics(*new Kinetics(setting)),
//   fieldForPcls(*new array4_double(
//     setting.grid().get_nxn(),
//     setting.grid().get_nyn(),
//     setting.grid().get_nzn(),
//     2*DFIELD_3or4)),
//   my_clock(0)
MIsolver::MIsolver(int argc, char **argv)
: setting(*new Setting(argc, argv)),
  miMoments(0),
  speciesMoms(0),
  EMf(0),
  speciesMoms(0),
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

  const int nxn = setting.grid().get_nxn();
  const int nyn = setting.grid().get_nyn();
  const int nzn = setting.grid().get_nzn();
  kinetics = new Kinetics(setting);
  speciesMoms = new SpeciesMoms(setting);
  miMoments = new MImoments(setting);
  EMf = new EMfields3D(setting, miMoments);
  fieldsForPcls = new array4_double(nxn,nyn,nzn,2*DFIELD_3or4);

  my_clock = new Timing(vct->get_rank());

  return 0;
}

void MIsolver::compute_moments()
{
  get_kinetics().compute_speciesMoms(speciesMoms);

  // miMoments are the moments needed to drive the implicit field
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
  miMoments->compute_from_speciesMoms(speciesMoms, Bx, By, Bz);
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
    EMf->calculateRhoHat(get_iMoments());
    // advance the E field
    EMf->calculateE(get_iMoments());
  }
  send_field_to_kinetic_solver();
}

void MIsolver::send_Bsmooth_to_kinetic_solver(
  arr3_double Bx_smooth,
  arr3_double By_smooth,
  arr3_double Bz_smooth)
{
  // send(EMf->fetch_Bx_smooth(),
  //      EMf->fetch_By_smooth(),
  //      EMf->fetch_Bz_smooth());
}

//! update Bfield (assuming Eth has already been calculated)
//  B^{n+1} = B^n - curl(Eth)
void MIsolver::advance_Bfield()
{
  timeTasks_set_main_task(TimeTasks::FIELDS);
  timeTasks_set_task(TimeTasks::BFIELD); // subtask
  // calculate the B field
  EMf->advanceB();

  if(I_am_field_solver())
  {
    // begin sending of B_smooth to kinetic solver
    send_Bsmooth_to_kinetic_solver();
  }
}

// this method should be a no-op on the cluster
void MIsolver::move_particles()
{
  if(I_am_kinetic_solver())
  {
    //[...receive field from fieldsolver...]
    kinetics->moveParticles();
  }
}

// === Section: general initialization_routines ===

/*! initialize Moments with initial configuration */
// todo: separate out this problem-specific code
//
void MIsolver::set_initial_conditions()
{
  string Case = col->getCase;
  int restart_status = col->getRestart_status();

  if (Case()=="Dipole")       initDipole(EMf, part);
  // cases prior to this point are responsible for their own restarts
  else if(restart_status)     init_from_restart(EMf, part);
  // cases that use rhon to set rhoc
  //else if (col->getCase()=="restart")   init_from_restart();
  else if (Case=="GEM")       initGEM(EMf, part);
  //else if (Case=="GEMnoPert") initGEMnoPert();
  else if (Case=="ForceFree") initForceFree(EMf, part);
  // cases that use rhoc to set rhon
#ifdef BATSRUS
  else if (Case=="BATSRUS")   initBATSRUS(EMf, part);
#endif
  //else if (Case=="RandomCase")initRandomField();
  else {
    eprintf("Case=%s in inputfile is not supported.", Case);
  }

  if(col->getRestart_status() == 0)
  {
    for (int i = 0; i < ns; i++)
    {
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

// sets:
//  * En, Bc, Bn, Bext, Btot,
//  * particles
void MIsolver::initialize()
{
  set_initial_conditions();
  // initialize total magnetic field
  EMf->update_total_B();
  // initialize moments from particles
  // (moments used to initialize particles are discarded).
  compute_moments();
  // write initial output before starting simulation
  initialize_output();
}

// === Section: specific initialization_routines ===

/*! initialize Moments with initial configuration */
void MIsolver::init_from_restart(EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  assert(col->getRestart_status()!=0)
  // READING FROM RESTART
  {
  #ifdef NO_HDF5
    eprintf("restart requires compiling with HDF5");
  #else
    read_moments_restart(&get_col(),vct,grid,&rhons,ns);

    // communicate species densities to ghost nodes
    for (int is = 0; is < ns; is++)
    {
      double ***moment0 = rhons.fetch_arr4()[is];
      communicateNode_P(nxn, nyn, nzn, moment0, vct);
    }

    if (col->getCase()=="Dipole") {
      ConstantChargePlanet(col->getL_square(),
        col->getx_center(),col->gety_center(),col->getz_center());
    }

    ConstantChargeOpenBC();
  #endif // NO_HDF5
  }

  for (int is = 0; is < ns; is++)
  {
    grid->interpN2C(rhocs, is, rhons);
    eprintf("need to read particles too. unimplemented.");
    //part[i].maxwellian(rhocs);
  }
}

/*! initialize Magnetic and Electric Field with initial configuration */
void MIsolver::init_from_restart(EMfields3D& EMf, Particles3Dcomm* part)
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  assert(col->getRestart_status());
  #ifdef NO_HDF5
  eprintf("restart requires compiling with HDF5");
  #else
  read_field_restart(&get_col(),vct,grid,Bxn,Byn,Bzn,Ex,Ey,Ez);
  {
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx, vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy, vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz, vct);

    // communicate E
    communicateNodeBC(nxn, nyn, nzn, Ex, col->bcEx, vct);
    communicateNodeBC(nxn, nyn, nzn, Ey, col->bcEy, vct);
    communicateNodeBC(nxn, nyn, nzn, Ez, col->bcEz, vct);
  }
  #endif // NO_HDF5
  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);

  // The ghost cells of Bxc are never used, so this is pointless. -eaj
  // communicate ghost
  //communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx, vct);
  //communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy, vct);
  //communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz, vct);
}

#ifdef BATSRUS
/*! initiliaze EM for GEM challange */
void MIsolver::initBATSRUS(EMfields3D& EMf, Particles3Dcomm* part)
{
  const Collective *col = &get_col();
  const Grid *grid = &get_grid();
  cout << "------------------------------------------" << endl;
  cout << "         Initialize from BATSRUS          " << endl;
  cout << "------------------------------------------" << endl;

  // populating these does nothing so has been removed.
  //
  //array4_double rhons(ns,nxn,nyn,nzn);
  //array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  // these are node-centered values.
  array3_double Exc(nxc,nyc,nzc);
  array3_double Eyc(nxc,nyc,nzc);
  array3_double Ezc(nxc,nyc,nzc);

  // loop over species and cell centers: fill in charge density
  //for (int is=0; is < ns; is++)
  //  for (int i=0; i < nxc; i++)
  //    for (int j=0; j < nyc; j++)
  //      for (int k=0; k < nzc; k++)
  //      {
  //        // WARNING getFluidRhoCenter contains "case" statment
  //        rhocs[is][i][j][k] = col->getFluidRhoCenter(i,j,k,is);
  //      }

  // loop over cell centers and fill in magnetic and electric fields
  for (int i=0; i < nxc; i++)
    for (int j=0; j < nyc; j++)
      for (int k=0; k < nzc; k++)
      {
        // This coupling was setting the electric field at
        // the nodes using Ohm's law applied to cell-centered
        // values, losing an order of accuracy.
        // I changed it to pass in cell-centered electric field
        // values and interpolate them to the nodes.
        // (The alternative would be to interpolate fluid
        // quantities to the nodes and apply Ohm's law
        // at the nodes.) -eaj
        //
        // WARNING getFluidRhoCenter contains "case" statment
        col->setFluidFieldsCenter(&Exc[i][j][k],&Eyc[i][j][k],&Ezc[i][j][k],
            &Bxc[i][j][k],&Byc[i][j][k],&Bzc[i][j][k],i,j,k);
      }

  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);
  grid->interpC2N(Ex,Exc);
  grid->interpC2N(Ey,Eyc);
  grid->interpC2N(Ez,Ezc);

  for (int i = 0; i < ns; i++)
  {
    //grid->interpC2N(rhons[is],rhocs[is]);
    part[i].MaxwellianFromFluid(col,i);
  }
}
#endif

/*! initiliaze EM for GEM challange */
void MIsolver::initGEM(EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  assert(col->getRestart_status()==0);
  // initialize
  {
    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
          Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXN(i, j, k) - Lx / 2;
          ypert = grid->getYN(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly );
          Byc[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXC(i, j, k) - Lx / 2;
          ypert = grid->getYC(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxc[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byc[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzc[i][j][k] = B0z;

        }
  }
  for (int i = 0; i < ns; i++)
  {
    grid->interpN2C(rhocs, is, rhons);
    part[i].maxwellian(rhocs);
  }
}

void c_Solver::initOriginalGEM()
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Bxn = fetch_Bxn();
  arr3_double Byn = fetch_Byn();
  arr3_double Bzn = fetch_Bzn();
  arr3_double Ex = fetch_Ex();
  arr3_double Ey = fetch_Ey();
  arr3_double Ez = fetch_Ez();
  const Grid *grid = &get_grid();
  const double Lx = get_col().getLx();
  const double Ly = get_col().getLy();
  const double Lz = get_col().getLz();
  const double B0x = get_col().getB0x();
  const double B0y = get_col().getB0y();
  const double B0z = get_col().getB0z();
  // initialize using perturbation localized in X
  {
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          const double yM = grid->getYN(i, j, k) - .5 * Ly;
          Bxn[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          Bxn[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byn[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          const double yM = grid->getYC(i, j, k) - .5 * Ly;
          Bxc[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXC(i, j, k) - .5 * Lx;
          Bxc[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byc[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzc[i][j][k] = B0z;
        }
  }
  for (int i = 0; i < ns; i++)
  {
    grid->interpN2C(rhocs, is, rhons);
    part[i].maxwellian(rhocs);
  }
}

void EMfields3D::set_Jext()
{
  // calculate the external current for reporting purposes
  //
  // set Bxn including Bx_ext
  for (int i=0; i < nxn; i++)
  for (int j=0; j < nyn; j++)
  for (int k=0; k < nzn; k++)
  {
    // We want a well-balanced sheme, where the equilibrium is an
    // exact solution of the discretized equations, so
    // Bx_ext should not be included in Bxn.
    // But we initially include it here as a means to
    // compute Jx_ext (to be used only for reporting purposes)
    // and then reset it.
    //
    Bxn[i][j][k] = B0x + fetch_Bx_ext()[i][j][k];
    Byn[i][j][k] = B0y + fetch_By_ext()[i][j][k];
    Bzn[i][j][k] = B0z + fetch_Bz_ext()[i][j][k];
  }
  //
  grid->interpN2C(Bxc,Bxn);
  grid->interpN2C(Byc,Byn);
  grid->interpN2C(Bzc,Bzn);
  //
  communicateCenterBC_P(nxc, nyc, nzc, Bxc, col->bcBx, vct);
  communicateCenterBC_P(nxc, nyc, nzc, Byc, col->bcBy, vct);
  communicateCenterBC_P(nxc, nyc, nzc, Bzc, col->bcBz, vct);
  //
  // initialize J_ext =c/4*pi curl(B) on nodes (current due to the dipole)
  //
  // the external current plays no part in the algorithm;
  // it is only for reporting the net current.
  assert(!Jx_ext);
  assert(!Jy_ext);
  assert(!Jz_ext);
  Jx_ext = new array3_double(nxn,nyn,nzn);
  Jy_ext = new array3_double(nxn,nyn,nzn);
  Jz_ext = new array3_double(nxn,nyn,nzn);
  grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc);
  scale(Jx_ext,tempXN,c/FourPI,nxn,nyn,nzn);
  scale(Jy_ext,tempYN,c/FourPI,nxn,nyn,nzn);
  scale(Jz_ext,tempZN,c/FourPI,nxn,nyn,nzn);
}

static void loopX(double *b, double z, double x, double y, double a,
  double zc, double xc, double yc, double m)
{
  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);
  //double Rho = r * sin(theta);
  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    eprintf("Err came back :%d", err);

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  //for debugging
  /*cout << "\n\nAt (" << x << "," << y << "," << z << "), the field is :" << endl;
    cout << "Bx: " << Bx << " T" << endl;
    cout << "By: " << By << " T" << endl;
    cout << "Bz: " << Bz << " T" << endl;
    cout << "BRho: " << BRho << " T" << endl;*/

  b[1] = Bx;
  b[2] = By;
  b[0] = Bz;
}

static void loopY(double *b, double y, double z, double x,
  double a, double yc, double zc, double xc, double m)
{
  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);
  //double Rho = r * sin(theta);
  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    eprintf("Err came back :%d", err);

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  //for debugging
  /*cout << "\n\nAt (" << x << "," << y << "," << z << "), the field is :" << endl;
    cout << "Bx: " << Bx << " T" << endl;
    cout << "By: " << By << " T" << endl;
    cout << "Bz: " << Bz << " T" << endl;
    cout << "BRho: " << BRho << " T" << endl;*/

  b[2] = Bx;
  b[0] = By;
  b[1] = Bz;
}

static void loopZ(double *b, double x, double y, double z,
  double a, double xc, double yc, double zc, double m)
{

  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);

  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    eprintf("Err came back :%d", err);

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  b[0] = Bx;
  b[1] = By;
  b[2] = Bz;
}

/*! Initialise a combination of magnetic dipoles */
//
// In this case, we run the initialization code even
// if we will overwrite the state data with restart data
// in order to make sure that Bext and Jext get initialized.
//
void MIsolver::initDipole(EMfields3D& EMf, Particles3Dcomm* part)
{
  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();
  arr3_double Bx_ext = EMf->fetch_Bx_ext();
  arr3_double By_ext = EMf->fetch_By_ext();
  arr3_double Bz_ext = EMf->fetch_Bz_ext();

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  const int B0x = col->getB0x();
  const int B0y = col->getB0y();
  const int B0z = col->getB0z();
  const int B1x = col->getB1x();
  const int B1y = col->getB1y();
  const int B1z = col->getB1z();

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  // initialize E
  for (int i=0; i < nxn; i++)
  for (int j=0; j < nyn; j++)
  for (int k=0; k < nzn; k++)
  {
    Ex[i][j][k] = ebc[0];
    Ey[i][j][k] = ebc[1];
    Ez[i][j][k] = ebc[2];
  }

  // Compute dipolar field B_ext
  //
  Bx_ext.setall(0.);
  By_ext.setall(0.);
  Bz_ext.setall(0.);
  // initialize external magnetic field
  for (int i=0; i < nxn; i++)
  for (int j=0; j < nyn; j++)
  for (int k=0; k < nzn; k++)
  {
    double blp[3];
    // Set coil diameter
    double a=delta;

    double xc=x_center;
    double yc=y_center;
    double zc=z_center;

    double x = grid->getXN(i,j,k);
    double y = grid->getYN(i,j,k);
    double z = grid->getZN(i,j,k);

    loopZ(blp, x, y, z, a, xc, yc, zc, B1z);
    Bx_ext[i][j][k]  += blp[0];
    By_ext[i][j][k]  += blp[1];
    Bz_ext[i][j][k]  += blp[2];
    loopX(blp, x, y, z, a, xc, yc, zc, B1x);
    Bx_ext[i][j][k] += blp[0];
    By_ext[i][j][k] += blp[1];
    Bz_ext[i][j][k] += blp[2];
    loopY(blp, x, y, z, a, xc, yc, zc, B1y);
    Bx_ext[i][j][k] += blp[0];
    By_ext[i][j][k] += blp[1];
    Bz_ext[i][j][k] += blp[2];
  }

  // this sets Bxn to include Bx_ext
  EMf->set_Jext();

  if (col->getRestart_status()==0)
  {
    // reset Bxn, excluding Bx_ext
    for (int i=0; i < nxn; i++)
    for (int j=0; j < nyn; j++)
    for (int k=0; k < nzn; k++)
    {
      Bxn[i][j][k] = B0x;
      Byn[i][j][k] = B0y;
      Bzn[i][j][k] = B0z;
    }
    // update Bxc
    {
      grid->interpN2C(Bxc,Bxn);
      grid->interpN2C(Byc,Byn);
      grid->interpN2C(Bzc,Bzn);
      //
      communicateCenterBC_P(nxc, nyc, nzc, Bxc, col->bcBx, vct);
      communicateCenterBC_P(nxc, nyc, nzc, Byc, col->bcBy, vct);
      communicateCenterBC_P(nxc, nyc, nzc, Bzc, col->bcBz, vct);
    }

    // set rhocs
    array4_double rhocs(ns,nxc,nyc,nzc);
    for (int is=0; is < ns; is++){
    for (int i=0; i < nxc; i++)
    for (int j=0; j < nyc; j++)
    for (int k=0; k < nzc; k++)
    {
      rhocs[is][i][j][k] = rhoINIT[is]/FourPI;
    }
    for (int is=0 ; is<ns; is++)
    {
      part[i].maxwellian(rhocs);
    }
  }
  else // assert(col->getRestart_status()!=0)
  {
    init_from_restart(EMf, part);  // use the fields from restart file
  }
}

// --- section: relegated_initializations ---
// initialization routines that should be relegated
// to a BloatedSolver that inherits from MIsolver
//#if 0 // other_initialization_routines

void MIsolver::initDoublePeriodicHarrisWithGaussianHumpPerturbation(
  EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Bxn = fetch_Bxn();
  arr3_double Byn = fetch_Byn();
  arr3_double Bzn = fetch_Bzn();
  arr3_double Ex = fetch_Ex();
  arr3_double Ey = fetch_Ey();
  arr3_double Ez = fetch_Ez();

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  const double Lx = get_col().getLx();
  const double Ly = get_col().getLy();
  const double Lz = get_col().getLz();
  const double B0x = get_col().getB0x();
  const double B0y = get_col().getB0y();
  const double B0z = get_col().getB0z();
  // perturbation localized in X
  const double pertX = 0.4;
  const double deltax = 8. * delta;
  const double deltay = 4. * delta;
  // initialize
  {
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              const double sech_yBd = 1. / cosh(yBd);
              const double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxn[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byn[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxn[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byn[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);

          // guide field
          Bzn[i][j][k] = B0z;
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx, vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy, vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz, vct);
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          Bxc[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxc[i][j][k] += 0.;
          Byc[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxc[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byc[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxc[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byc[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);
          // guide field
          Bzc[i][j][k] = B0z;
        }
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx, vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy, vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz, vct);
    for (int is = 0; is < ns; is++)
    {
      grid->interpN2C(rhocs, is, rhons);
      part[i].maxwellian(rhocs);
    }
  }
}


/*! initialize GEM challenge with no Perturbation with dipole-like tail topology */
void MIsolver::initGEMDipoleLikeTailNoPert(
  EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Bxn = fetch_Bxn();
  arr3_double Byn = fetch_Byn();
  arr3_double Bzn = fetch_Bzn();
  arr3_double Ex = fetch_Ex();
  arr3_double Ey = fetch_Ey();
  arr3_double Ez = fetch_Ez();

  const Grid *grid = &get_grid();
  const double Lx = get_col().getLx();
  const double Ly = get_col().getLy();
  const double Lz = get_col().getLz();
  const double B0x = get_col().getB0x();
  const double B0y = get_col().getB0y();
  const double B0z = get_col().getB0z();
  const Grid *grid = &get_grid();
  // parameters controling the field topology
  // e.g., x1=Lx/5,x2=Lx/4 give 'separated' fields, x1=Lx/4,x2=Lx/3 give 'reconnected' topology

  double x1 = Lx / 6.0;         // minimal position of the gaussian peak 
  double x2 = Lx / 4.0;         // maximal position of the gaussian peak (the one closer to the center)
  double sigma = Lx / 15;       // base sigma of the gaussian - later it changes with the grid
  double stretch_curve = 2.0;   // stretch the sin^2 function over the x dimension - also can regulate the number of 'knots/reconnecitons points' if less than 1
  double skew_parameter = 0.50; // skew of the shape of the gaussian
  double pi = 3.1415927;
  double r1, r2, delta_x1x2;

  assert(col->getRestart_status()==0);
  // initialize
  {
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYN(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYN(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYC(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYC(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
    {
      grid->interpN2C(rhocs, is, rhons);
      part[i].maxwellian(rhocs);
    }
  }
}

/*! initialize GEM challenge with no Perturbation */
void MIsolver::initGEMnoPert(EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Bxn = fetch_Bxn();
  arr3_double Byn = fetch_Byn();
  arr3_double Bzn = fetch_Bzn();
  arr3_double Ex = fetch_Ex();
  arr3_double Ey = fetch_Ey();
  arr3_double Ez = fetch_Ez();

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  assert(col->getRestart_status()==0);
  // initialize
  {
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
    {
      grid->interpN2C(rhocs, is, rhons);
      part[i].maxwellian(rhocs);
    }
  }
}

// new init, random problem
void MIsolver::initRandomField(EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  double **modes_seed = newArr2(double, 7, 7);
  assert(col->getRestart_status()!=0)
  // initialize
  {
    if (get_vct().getCartesian_rank() ==0){
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i=0; i < ns; i++){
	cout << "rho species " << i <<" = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    double kx;
    double ky;
        
    /*       stringstream num_proc;
	     num_proc << vct->getCartesian_rank() ;
	     string cqsat = SaveDirName + "/RandomNumbers" + num_proc.str() + ".txt";
        ofstream my_file(cqsat.c_str(), fstream::binary);
	for (int m=-3; m < 4; m++)
            for (int n=-3; n < 4; n++){
            modes_seed[m+3][n+3] = rand() / (double) RAND_MAX;
            my_file <<"modes_seed["<< m+3<<"][" << "\t" << n+3 << "] = " << modes_seed[m+3][n+3] << endl;
            }
              my_file.close();
    */
    modes_seed[0][0] = 0.532767;
    modes_seed[0][1] = 0.218959;
    modes_seed[0][2] = 0.0470446;
    modes_seed[0][3] = 0.678865;
    modes_seed[0][4] = 0.679296;
    modes_seed[0][5] = 0.934693;
    modes_seed[0][6] = 0.383502;
    modes_seed[1][0] = 0.519416;
    modes_seed[1][1] = 0.830965;
    modes_seed[1][2] = 0.0345721;
    modes_seed[1][3] = 0.0534616;
    modes_seed[1][4] = 0.5297;
    modes_seed[1][5] = 0.671149;
    modes_seed[1][6] = 0.00769819;
    modes_seed[2][0] = 0.383416;
    modes_seed[2][1] = 0.0668422;
    modes_seed[2][2] = 0.417486;
    modes_seed[2][3] = 0.686773;
    modes_seed[2][4] = 0.588977;
    modes_seed[2][5] = 0.930436;
    modes_seed[2][6] = 0.846167;
    modes_seed[3][0] = 0.526929;
    modes_seed[3][1] = 0.0919649;
    modes_seed[3][2] = 0.653919;
    modes_seed[3][3] = 0.415999;
    modes_seed[3][4] = 0.701191;
    modes_seed[3][5] = 0.910321;
    modes_seed[3][6] = 0.762198;
    modes_seed[4][0] = 0.262453;
    modes_seed[4][1] = 0.0474645;
    modes_seed[4][2] = 0.736082;
    modes_seed[4][3] = 0.328234;
    modes_seed[4][4] = 0.632639;
    modes_seed[4][5] = 0.75641;
    modes_seed[4][6] = 0.991037;
    modes_seed[5][0] = 0.365339;
    modes_seed[5][1] = 0.247039;
    modes_seed[5][2] = 0.98255;
    modes_seed[5][3] = 0.72266;
    modes_seed[5][4] = 0.753356;
    modes_seed[5][5] = 0.651519;
    modes_seed[5][6] = 0.0726859;
    modes_seed[6][0] = 0.631635;
    modes_seed[6][1] = 0.884707;
    modes_seed[6][2] = 0.27271;
    modes_seed[6][3] = 0.436411;
    modes_seed[6][4] = 0.766495;
    modes_seed[6][5] = 0.477732;
    modes_seed[6][6] = 0.237774;

    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++){
	  // initialize the density for species
	  for (int is=0; is < ns; is++){
	    rhons[is][i][j][k] = rhoINIT[is]/FourPI;
	  }
	  // electric field
	  Ex[i][j][k] =  0.0;
	  Ey[i][j][k] =  0.0;
	  Ez[i][j][k] =  0.0;
	  // Magnetic field
	  Bxn[i][j][k] =  0.0;
	  Byn[i][j][k] =  0.0;
	  Bzn[i][j][k] =  B0z;
	  for (int m=-3; m < 4; m++)
	    for (int n=-3; n < 4; n++){

	      kx=2.0*M_PI*m/Lx;
	      ky=2.0*M_PI*n/Ly;
	      Bxn[i][j][k] += -B0x*ky*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	      Byn[i][j][k] += B0x*kx*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	      // Bzn[i][j][k] += B0x*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	    }
	}
	  // communicate ghost
	  communicateNodeBC(nxn, nyn, nzn, Bxn, 1, 1, 2, 2, 1, 1, vct);
	  communicateNodeBC(nxn, nyn, nzn, Byn, 1, 1, 1, 1, 1, 1, vct);
	  communicateNodeBC(nxn, nyn, nzn, Bzn, 1, 1, 2, 2, 1, 1, vct);
	  // initialize B on centers
	  grid->interpN2C(Bxc, Bxn);
	  grid->interpN2C(Byc, Byn);
	  grid->interpN2C(Bzc, Bzn);
	  // communicate ghost
	  communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct);
	  communicateCenterBC(nxc, nyc, nzc, Byc, 1, 1, 1, 1, 1, 1, vct);
	  communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct);
	  for (int is=0 ; is<ns; is++)
            grid->interpN2C(rhocs,is,rhons);
  }
  delArr2(modes_seed, 7);
}


/*! Init Force Free (JxB=0) */
void MIsolver::initForceFree(EMfields3D& EMf, Particles3Dcomm* part)
{
  array4_double rhons(ns,nxn,nyn,nzn);
  array4_double rhocs(ns,nxc,nyc,nzc);

  arr3_double Ex = EMf->fetch_Ex();
  arr3_double Ey = EMf->fetch_Ey();
  arr3_double Ez = EMf->fetch_Ez();
  arr3_double Bxc = EMf->fetch_Bxc();
  arr3_double Byc = EMf->fetch_Byc();
  arr3_double Bzc = EMf->fetch_Bzc();
  arr3_double Bxn = EMf->fetch_Bxn();
  arr3_double Byn = EMf->fetch_Byn();
  arr3_double Bzn = EMf->fetch_Bzn();

  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  assert (col->getRestart_status() == 0);
  // initialize
  {
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------" << endl;
      cout << "Initialize Force Free with Perturbation" << endl;
      cout << "----------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
      }
      cout << "Smoothing Factor = " << Smooth << endl;
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          Bxn[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXN(i, j, k) / Lx) * sin(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          Byn[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXN(i, j, k) / Lx) * cos(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzn[i][j][k] = B0z / cosh((grid->getYN(i, j, k) - Ly / 2) / delta);
        }
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the perturbation
          Bxc[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXC(i, j, k) / Lx) * sin(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          Byc[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXC(i, j, k) / Lx) * cos(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzc[i][j][k] = B0z / cosh((grid->getYC(i, j, k) - Ly / 2) / delta);
        }

    for (int is = 0; is < ns; is++)
    {
      grid->interpN2C(rhocs, is, rhons);
      part[i].force_free(is, rhocs);
    }
  }
}
//#endif // other_initialization_routines

// --- end of relegated_initializations ---
// === end of specific initialization_routines ===

// === section: output methods ===

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
// C    En          = advance(dt, En, Bn, Btot, Jhat, rhons)
// C    Bc          = advance(dt, Bc, En)
// C    Bn          = interpolate_to_nodes(Bc)
// C    Btot        = Bn + Bext
// C B  B_smooth    = smooth(Btot)
// C B  E_smooth    = smooth(En)
//   B  BEaos       = soa2aos(B_smooth,E_smooth)
//   B  particles   = advance(dt, particles, BEaos)
// C B  speciesMoms    = sumMoments(particles)
// C    Jhat_coarse = compute_Jhat(speciesMoms, B_smooth)
// C    Jhat        = smooth(compute_Jhat(speciesMoms, B_smooth)
// C    rhons       = smooth(speciesMoms.rhons)
//
// remarks:
//  
// * we could transfer Jhat_coarse rather than speciesMoms
// * a more correct way to compute Jhat is:
//   C    Jhat        = compute_Jhat(smooth(SpeciesMoms), Btot)
// * Btot is used to compute MUdot in En.advance().
// * Btot or B_smooth is used in Jhat.smooth() to compute PIdot.
// * Jhat could be transferred from B to C rather than SpeciesMoms.
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
//   speciesMoms
//   particles
// FieldSolver:
//   EMfields3D:
//     En, Bc, Bn, Btot
//   B_smooth
//   E_smooth
//   speciesMoms
//   Jhat_coarse
//   Jhat
//   rhon
// KineticSolver:
//   B_smooth
//   E_smooth
//   BEaos
//   particles
//   speciesMoms
//
void MIsolver::run()
{
  timeTasks.resetCycle();
  initialize();
  // shouldn't we call this?
  //WriteOutput(FirstCycle()-1);
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
