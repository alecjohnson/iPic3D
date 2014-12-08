#include "Moments.h"
#include "Alloc.h"

void Moments10::set_to_zero()
{
  arr.setall(0);
  //#pragma omp parallel for collapse(4)
  //for (register int i = 0; i < nx; i++)
  //for (register int j = 0; j < ny; j++)
  //for (register int k = 0; k < nz; k++)
  //for (register int m = 0; m < 10; m++)
  //{
  //  arr[i][j][k][m] = 0.0;
  //}
}

SpeciesMoms::~SpeciesMoms()
{
  for(int i=0;i<sizeMomentsArray;i++) { delete moments10Array[i]; }
  delete [] moments10Array;
}

ReportedMoments::ReportedMoments(const Setting& setting_)
: setting(setting_),
  nxc(setting.grid().get_nxc()),
  nyc(setting.grid().get_nyc()),
  nzc(setting.grid().get_nzc()),
  nxn(setting.grid().get_nxn()),
  nyn(setting.grid().get_nyn()),
  nzn(setting.grid().get_nzn()),
  // array allocation: nodes
  //
  Jx   (nxn, nyn, nzn),
  Jy   (nxn, nyn, nzn),
  Jz   (nxn, nyn, nzn),
  //Jxh  (nxn, nyn, nzn),
  //Jyh  (nxn, nyn, nzn),
  //Jzh  (nxn, nyn, nzn),
  //
{
}

SpeciesMoms::SpeciesMoms(const Setting& setting_, arr4_double rhons_)
: setting(setting_),
  grid(&setting.get_grid()),
  col(&setting.get_col()),
  nxc(grid->get_nxc()),
  nyc(grid->get_nyc()),
  nzc(grid->get_nzc()),
  nxn(grid->get_nxn()),
  nyn(grid->get_nyn()),
  nzn(grid->get_nzn()),
  invVOL(grid->getInvVOL()),
  ns(col->getNs()),
  //
  // species-specific quantities
  //
  //rhons (ns, nxn, nyn, nzn),
  //rhocs (ns, nxc, nyc, nzc),
  // external owned memory
  rhons(rhons_),
  // internally owned memory
  Jxs   (ns, nxn, nyn, nzn),
  Jys   (ns, nxn, nyn, nzn),
  Jzs   (ns, nxn, nyn, nzn),
  pXXsn (ns, nxn, nyn, nzn),
  pXYsn (ns, nxn, nyn, nzn),
  pXZsn (ns, nxn, nyn, nzn),
  pYYsn (ns, nxn, nyn, nzn),
  pYZsn (ns, nxn, nyn, nzn),
  pZZsn (ns, nxn, nyn, nzn),

  // array allocation: central points 
  //
  //PHI  (nxc, nyc, nzc),
  //Bxc  (nxc, nyc, nzc),
  //Byc  (nxc, nyc, nzc),
  //Bzc  (nxc, nyc, nzc),
  //rhoc (nxc, nyc, nzc),

  // temporary arrays to calculate hatted moments
  //
  tempXC (nxc, nyc, nzc),
  tempYC (nxc, nyc, nzc),
  tempZC (nxc, nyc, nzc),
  //
  tempXN (nxn, nyn, nzn),
  tempYN (nxn, nyn, nzn),
  tempZN (nxn, nyn, nzn)
{
  if(Parameters::get_VECTORIZE_MOMENTS())
  {
    // In this case particles are sorted
    // and there is no need for each thread
    // to sum moments in a separate array.
    sizeMomentsArray = 1;
  }
  else
  {
    sizeMomentsArray = omp_get_max_threads();
  }
  moments10Array = new Moments10*[sizeMomentsArray];
  for(int i=0;i<sizeMomentsArray;i++)
  {
    moments10Array[i] = new Moments10(nxn,nyn,nzn);
  }
}

// === section: methods_to_compute_implicit_moments ===

/*! Calculate PI dot (vectX, vectY, vectZ) */
void SpeciesMoms::PIdot(
  arr3_double PIdotX,
  arr3_double PIdotY,
  arr3_double PIdotZ,
  const_arr3_double vectX,
  const_arr3_double vectY,
  const_arr3_double vectZ,
  const_arr3_double Bx,
  const_arr3_double By,
  const_arr3_double Bz,
  double qom_)
{
  const Grid *grid = &get_grid();
  double beta, edotb, omcx, omcy, omcz, denom;
  beta = .5 * qom_ * dt / c;
  for (int i = 1; i < nxn - 1; i++)
  for (int j = 1; j < nyn - 1; j++)
  for (int k = 1; k < nzn - 1; k++)
  {
    omcx = beta * Bx[i][j][k];
    omcy = beta * By[i][j][k];
    omcz = beta * Bz[i][j][k];
    edotb = vectX.get(i,j,k) * omcx + vectY.get(i,j,k) * omcy + vectZ.get(i,j,k) * omcz;
    denom = 1 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
    PIdotX.fetch(i,j,k) += (vectX.get(i,j,k) + (vectY.get(i,j,k) * omcz - vectZ.get(i,j,k) * omcy + edotb * omcx)) * denom;
    PIdotY.fetch(i,j,k) += (vectY.get(i,j,k) + (vectZ.get(i,j,k) * omcx - vectX.get(i,j,k) * omcz + edotb * omcy)) * denom;
    PIdotZ.fetch(i,j,k) += (vectZ.get(i,j,k) + (vectX.get(i,j,k) * omcy - vectY.get(i,j,k) * omcx + edotb * omcz)) * denom;
  }
}

/*! Calculate Jx hat, Jy hat, Jz hat */
void SpeciesMoms::calculateJhat(
  arr3_double Jxh,
  arr3_double Jyh,
  arr3_double Jzh,
  const_arr3_double Bx,
  const_arr3_double By,
  const_arr3_double Bz)
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  // temporary arrays to compute hatted moments
  //
  array3_double tempXC (nxc, nyc, nzc);
  array3_double tempYC (nxc, nyc, nzc);
  array3_double tempZC (nxc, nyc, nzc);
  array3_double tempXN (nxn, nyn, nzn);
  array3_double tempYN (nxn, nyn, nzn);
  array3_double tempZN (nxn, nyn, nzn);

  Jxh.setall(0.);
  Jyh.setall(0.);
  Jzh.setall(0.);

  for (int is = 0; is < ns; is++) {
    grid->divSymmTensorN2C(tempXC, tempYC, tempZC,
      pXXsn, pXYsn, pXZsn, pYYsn, pYZsn, pZZsn, is);

    scale(tempXC, -dt / 2.0, nxc, nyc, nzc);
    scale(tempYC, -dt / 2.0, nxc, nyc, nzc);
    scale(tempZC, -dt / 2.0, nxc, nyc, nzc);
    // communicate before interpolating
    communicateCenterBC_P(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
    communicateCenterBC_P(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
    communicateCenterBC_P(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);

    grid->interpC2N(tempXN, tempXC);
    grid->interpC2N(tempYN, tempYC);
    grid->interpC2N(tempZN, tempZC);
    sum(tempXN, Jxs, nxn, nyn, nzn, is);
    sum(tempYN, Jys, nxn, nyn, nzn, is);
    sum(tempZN, Jzs, nxn, nyn, nzn, is);

    if(Parameters::use_perfect_smoothing())
    {
      // smooth before applying magnetic field
      // (requires smoothing num_species times as many quantities
      // as if we wait to smooth until after computing Jhat).
      for(int i=0; i<Parameters::get_num_smoothings(); i++);
      {
        get_grid().smooth(tempXN, 1);
        get_grid().smooth(tempYN, 1);
        get_grid().smooth(tempZN, 1);
      }
    }
    PIdot(Jxh, Jyh, Jzh, tempXN, tempYN, tempZN,
      Bx, By, Bz, qom[is])
  }
  // This original way smooths the current only once; this
  // could explain why it becomes necessary to retain
  // smoothing in the electric field.
  if(Parameters::use_original_smoothing())
  {
    get_grid().smooth(Jxh, 1);
    get_grid().smooth(Jyh, 1);
    get_grid().smooth(Jzh, 1);
  }
  else if(Parameters::use_correct_smoothing()
      && !Parameters::use_perfect_smoothing())
  {
    for(int i=0; i<Parameters::get_num_smoothings(); i++);
    {
      get_grid().smooth(Jxh, 1);
      get_grid().smooth(Jyh, 1);
      get_grid().smooth(Jzh, 1);
    }
  }
}

// === section: methods_to_accumulate_moments ===
//
// These methods do not perform any communication,
// so boundary nodes contain only the portion that
// this process contributes to.
//
void SpeciesMoms::accumulateMoments(const int is, const Particles3Dcomm& pcls)
{
  assert_eq(is, pcls.get_species_num());

  pcls.pad_capacities();

  if(Parameters::get_VECTORIZE_MOMENTS())
  {
    setZeroSpeciesMoms(is);
    sumMoments_vec(pcls);
    // these assume that particles are sorted by mesh cell
    //switch(Parameters::get_MOMENTS_TYPE())
    //{
    //  case Parameters::SoA:
    //    // since particles are sorted,
    //    // we can vectorize interpolation of particles to grid
    //    convertParticlesToSoA();
    //    pcls.sortParticles();
    //    sumMoments_vectorized(pcls);
    //    break;
    //  case Parameters::AoS:
    //    convertParticlesToAoS();
    //    sortParticles();
    //    sumMoments_vectorized_AoS(pcls);
    //    break;
    //  default:
    //    unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    //}
  }
  else
  {
    setZeroSpeciesMoms(is);
    if(Parameters::get_SORTING_PARTICLES())
      pcls.sortParticles();
    switch(Parameters::get_MOMENTS_TYPE())
    {
      case Parameters::SoA:
        pcls.convertParticlesToSoA();
        sumMoments(pcls);
        //sumMomentsOld(pcls);
        break;
      case Parameters::AoS:
        pcls.convertParticlesToAoS();
        sumMoments_AoS(part);
        break;
      case Parameters::AoSintr:
        pcls.convertParticlesToAoS();
        sumMoments_AoS_intr(part);
        break;
      default:
        unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    }
  }
}

// This was Particles3Dcomm::interpP2G()
void SpeciesMoms::sumMomentsOld(const Particles3Dcomm& pcls)
{
  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / grid->getDX();
  const double inv_dy = 1.0 / grid->getDY();
  const double inv_dz = 1.0 / grid->getDZ();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  double const*const x = pcls.getXall();
  double const*const y = pcls.getYall();
  double const*const z = pcls.getZall();
  double const*const u = pcls.getUall();
  double const*const v = pcls.getVall();
  double const*const w = pcls.getWall();
  double const*const q = pcls.getQall();
  //
  const int is = pcls.get_species_num();

  const int nop = pcls.getNOP();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP
  TimeTasks timeTasksAcc;
  #pragma omp parallel private(timeTasks)
  {
    int thread_num = omp_get_thread_num();
    timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION);
    Moments10& speciesMoments10 = fetch_moments10Array(thread_num);
    speciesMoments10.set_to_zero();
    arr4_double moments = speciesMoments10.fetch_arr();
    // The following loop is expensive, so it is wise to assume that the
    // compiler is stupid.  Therefore we should on the one hand
    // expand things out and on the other hand avoid repeating computations.
    #pragma omp for
    for (int i = 0; i < nop; i++)
    {
      // compute the quadratic moments of velocity
      //
      const double ui=u[i];
      const double vi=v[i];
      const double wi=w[i];
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[10];
      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;

      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      const double xi0   = x[i] - grid->getXN(ix-1);
      const double eta0  = y[i] - grid->getYN(iy-1);
      const double zeta0 = z[i] - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - x[i];
      const double eta1  = grid->getYN(iy) - y[i];
      const double zeta1 = grid->getZN(iz) - z[i];
      const double qi = q[i];
      const double weight000 = qi * xi0 * eta0 * zeta0 * invVOL;
      const double weight001 = qi * xi0 * eta0 * zeta1 * invVOL;
      const double weight010 = qi * xi0 * eta1 * zeta0 * invVOL;
      const double weight011 = qi * xi0 * eta1 * zeta1 * invVOL;
      const double weight100 = qi * xi1 * eta0 * zeta0 * invVOL;
      const double weight101 = qi * xi1 * eta0 * zeta1 * invVOL;
      const double weight110 = qi * xi1 * eta1 * zeta0 * invVOL;
      const double weight111 = qi * xi1 * eta1 * zeta1 * invVOL;
      double weights[8];
      weights[0] = weight000;
      weights[1] = weight001;
      weights[2] = weight010;
      weights[3] = weight011;
      weights[4] = weight100;
      weights[5] = weight101;
      weights[6] = weight110;
      weights[7] = weight111;

      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        momentsArray[0] = moments[ix  ][iy  ][iz  ]; // moments000 
        momentsArray[1] = moments[ix  ][iy  ][iz-1]; // moments001 
        momentsArray[2] = moments[ix  ][iy-1][iz  ]; // moments010 
        momentsArray[3] = moments[ix  ][iy-1][iz-1]; // moments011 
        momentsArray[4] = moments[ix-1][iy  ][iz  ]; // moments100 
        momentsArray[5] = moments[ix-1][iy  ][iz-1]; // moments101 
        momentsArray[6] = moments[ix-1][iy-1][iz  ]; // moments110 
        momentsArray[7] = moments[ix-1][iy-1][iz-1]; // moments111 

        for(int m=0; m<10; m++)
        for(int c=0; c<8; c++)
        {
          momentsArray[c][m] += velmoments[m]*weights[c];
        }
      }
    }
    timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

    // reduction
    timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

    // reduce arrays
    {
      #pragma omp critical (reduceMoment0)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { rhons[is][i][j][k] += invVOL*moments[i][j][k][0]; }}
      #pragma omp critical (reduceMoment1)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1]; }}
      #pragma omp critical (reduceMoment2)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jys  [is][i][j][k] += invVOL*moments[i][j][k][2]; }}
      #pragma omp critical (reduceMoment3)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3]; }}
      #pragma omp critical (reduceMoment4)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4]; }}
      #pragma omp critical (reduceMoment5)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5]; }}
      #pragma omp critical (reduceMoment6)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6]; }}
      #pragma omp critical (reduceMoment7)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7]; }}
      #pragma omp critical (reduceMoment8)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8]; }}
      #pragma omp critical (reduceMoment9)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9]; }}
    }
    timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
    #pragma omp critical
    timeTasksAcc += timeTasks;
  }
  // reset timeTasks to be its average value for all threads
  timeTasksAcc /= omp_get_max_threads();
  timeTasks = timeTasksAcc;
  //communicateGhostP2G(is);
}
//
// Create a vectorized version of the moment accumulator as follows.
//
// A. Moment accumulation
//
// Let N:=sizeof(vector_unit)/sizeof(double).
// Case 1: Assuming AoS particle layout and using intrinsics vectorization:
//   Process P:=N/4 particles at a time:
//   1. gather position coordinates from P particles and
//      generate Px8 array of weights and P cell indices.
//   2. for each particle, add 10x8 array of moment-weight
//      products to appropriate cell accumulator.
//   Each cell now has a 10x8 array of node-destined moments.
//   (See sumMoments_AoS_intr().)
// Case 2: Assuming SoA particle layout and using trivial vectorization:
//   Process N particles at a time:
//   1. for pcl=1:N: (3) positions -> (8) weights, cell_index
//   2. For each of 10 moments:
//      a. for pcl=1:N: (<=2 of 3) charge velocities -> (1) moment
//      b. for pcl=1:N: (1) moment, (8) weights -> (8) node-destined moments
//      c. transpose 8xN array of node-destined moments to Nx8 array 
//      d. foreach pcl: add node-distined moments to cell of cell_index
//   Each cell now has a 10x8 array of node-destined moments.
//
//   If particles are sorted by mesh cell, then all moments are destined
//   for the same node; in this case, we can simply accumulate an 8xN
//   array of node-destined moments in each mesh cell and at the end
//   gather these moments at the nodes; to help performance and
//   code reuse, we will in each cell first transpose the 8xN array
//   of node-destined moments to an Nx8 array.
//
// B: Moment reduction
//
//   Gather the information from cells to nodes:
//   3. [foreach cell transpose node-destined moments:
//      10x8 -> 8x10, or rather (8+2)x8 -> 8x8 + 8x2]
//   4. at each node gather moments from cells.
//   5. [transpose moments at nodes if step 3 was done.]
//
//   We will likely omit steps 3 and 5; they could help to optimize,
//   but even without these steps, step 4 is not expected to dominate
//   unless the number of particles per mesh cell is small.
//
void SpeciesMoms::sumMoments_vec(const Particles3Dcomm& pcls)
{
  // Process N:=sizeof(vector_unit)/sizeof(double) particles at a time:
  const int Np=8;
  const int num_threads = omp_get_max_threads();

  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / grid->getDX();
  const double inv_dy = 1.0 / grid->getDY();
  const double inv_dz = 1.0 / grid->getDZ();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const int nxc = grid->getNXC();
  const int nyc = grid->getNYC();
  const int nzc = grid->getNZC();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  //
  // allocate memory to accumulate moments in each cell that
  // are destined for the nodes. we need to include ghost
  // cells even though nothing should go in them so that the
  // code that gathers moments will work.
  const int ncells=nxc*nyc*nzc;
  #pragma omp single
  double ****node_destined_moms = newArray4<double>(num_threads,10,ncells,8);
  // if we want to vectorize the reduction, we would need
  // to dimension this as:
  //double ****node_destined_moms
  //  = newArr4<double>(num_threads,ncells,10,8);
  // prices that we would pay:
  // - need to transpose prior to gathering moments.
  // - need to transpose moments again before using them in field solver.
  // limited benefits forseen:
  // - vectorized transpose requires handling the last two
  //   moments separately with basically scalar performance.
  //
  // indexing convention:
  //   it = thread index,    0:num_threads-1
  //   im = moment index,    0:10-1
  //   ic = cell index,      0:ncells-1
  //   id = direction index, 0:8-1
  //   ip = particle index,  0:Np-1, Np=8
  //
  {
      const int is = pcls.get_species_num();

      const vector_SpeciesParticle& pcl_list = pcls.get_pcl_list();

      const int nop = pcls.getNOP();

      int thread_num = omp_get_thread_num();
      if(!thread_num) { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
      //
      double ***node_destined_moms_arr = node_destined_moms[thread_num];
      // clear moments array
      {
	double *arrptr = &node_destined_moms_arr[0][0][0];
	const int numel = ncells*10*8;
	#pragma simd // this should vectorize
        for(int i=0; i<numel;i++) arrptr[i]=0;
      }
      #pragma omp for
      for (int p0 = 0; p0 < nop; p0+=Np)
      {
        const double *u,*v,*w,*q;
        const double *x,*y,*z;
        double storage_for_AoS_case[8][Np] ALLOC_ALIGNED;
        // if SoA
        if(pcls.get_particleType() == ParticleType::SoA)
        {
          u = &(pcls.getUall()[p0]);
          v = &(pcls.getVall()[p0]);
          w = &(pcls.getWall()[p0]);
          q = &(pcls.getQall()[p0]);
          x = &(pcls.getXall()[p0]);
          y = &(pcls.getYall()[p0]);
          z = &(pcls.getZall()[p0]);
        }
        // else convert this block to AoS
        else if(pcls.get_particleType() == ParticleType::AoS)
        {
          // use fast matrix transpose to generate block to process
          assert_eq(Np,8);
          const double(*in)[8]=reinterpret_cast<const double(*)[8]>(&pcl_list[p0]);
          double(*out)[8] = reinterpret_cast<double(*)[8]>(storage_for_AoS_case[0]);
          transpose_8x8_double(storage_for_AoS_case, in);
          u = storage_for_AoS_case[0];
          v = storage_for_AoS_case[1];
          w = storage_for_AoS_case[2];
          q = storage_for_AoS_case[3];
          x = storage_for_AoS_case[4];
          y = storage_for_AoS_case[5];
          z = storage_for_AoS_case[6];
        }
        else
        {
          unsupported_value_error(pcls.get_particleType());
        }
        ASSUME_ALIGNED(u);
        ASSUME_ALIGNED(v);
        ASSUME_ALIGNED(w);
        ASSUME_ALIGNED(q);
        ASSUME_ALIGNED(x);
        ASSUME_ALIGNED(y);
        ASSUME_ALIGNED(z);

        // 1. for pcl=1:Np: (3) positions -> (8) weights, cell_index
        double weights[8][Np] ALLOC_ALIGNED;
        int cell_index[Np]; // one-dimensional index of underlying array
        // will the compiler be smart enough to expand and vectorize this?
        #pragma simd
        for(int ip=0;ip<Np;ip++)
        {
          int cx,cy,cz;
          grid->get_cell_and_weights(
            x[ip], y[ip], z[ip],
            cx,cy,cz,
            weights[0][ip],
            weights[1][ip],
            weights[2][ip],
            weights[3][ip],
            weights[4][ip],
            weights[5][ip],
            weights[6][ip],
            weights[7][ip]);
          cell_index[ip] = grid->get_cell_index(cx,cy,cz);
        }

        // 2. For each of 10 moments:
        //    a. for pcl=1:Np: (<=2 of 3) charge velocities -> (1) moment
        double velmoments[10][Np];
        #pragma simd
        for(int ip=0;ip<Np;ip++)
        {
          const double ui=u[ip];
          const double vi=v[ip];
          const double wi=w[ip];
          velmoments[0][ip] = 1.;
          velmoments[1][ip] = ui;
          velmoments[2][ip] = vi;
          velmoments[3][ip] = wi;
          velmoments[4][ip] = ui*ui;
          velmoments[5][ip] = ui*vi;
          velmoments[6][ip] = ui*wi;
          velmoments[7][ip] = vi*vi;
          velmoments[8][ip] = vi*wi;
          velmoments[9][ip] = wi*wi;
        }
        // double node_moments[10][8][Np];
        for(int im=0;im<10;im++)
        {
          double **node_moms_arr = node_destined_moms_arr[im];
          //  b. for pcl=1:Np: (1)moment, (8)weights -> (8)node-destined moments
          double node_moments[8][Np];
          for(int id=0;id<8;id++)
          #pragma simd
          for(int ip=0;ip<Np;ip++)
          {
            node_moments[id][ip]=weights[id][ip]*velmoments[im][ip];
          }
          //  c. transpose 8xN array of node-destined moments to Nx8 array 
          assert_eq(Np,8);
          transpose_8x8_double(node_moments);
          //  d. foreach pcl: add node-distined moments to cell of cell_index
          for(int ip=0;ip<Np;ip++)
          {
            double* cell_node_moms = node_moms_arr[cell_index[ip]];
            double* node_moms = node_moments[ip];
            #pragma simd
            for(int id=0;id<8;id++)
            {
              cell_node_moms[id] += node_moms[id];
            }
          }
        }
      }
      // Each cell now has a 10x8 array of node-destined moments.
      // Compute the quadratic moments of velocity
      //
      if(!thread_num) timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

      // reduction
      if(!thread_num) timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

      arr_fetch3(double) out_moments[10] =
      {
        rhons[is],
        Jxs  [is],
        Jys  [is],
        Jzs  [is],
        pXXsn[is],
        pXYsn[is],
        pXZsn[is],
        pYYsn[is],
        pYZsn[is],
        pZZsn[is]
      };

      // reduce moments in parallel.
      //
      // To vectorize this reduction, we
      // would need to change the indices of
      // the node-destined moments array from
      // node_destined_moms[it][im][ic][id] to
      // node_destined_moms[it][ic][im][id],
      // transpose the last two indices of
      // node_destined_moms and then gather.
      // the field solver would then need to
      // transpose the moments at the nodes.
      //
      if(true)
      {
        // this gathers node data from neighboring cells.
        // note that ghost nodes should not get any moments,
        // but we need them to exist for this code to work.
        // could add a check to verify that ghost node
        // moments are zero.
	//
        #pragma omp for
        for(int it=0;it<num_threads;it++)
        for(int im=0;im<10;im++)
        for(int ix=1;ix<nxc;ix++)
        for(int iy=1;iy<nyc;iy++)
        for(int iz=1;iz<nzc;iz++)
        {
          const int cx=ix-1;
          const int cy=iy-1;
          const int cz=iz-1;
          // gather data from neighboring cells intended for this node
          out_moments[im][ix][iy][iz] += invVOL*(
            node_destined_moms[it][im][grid->get_cell_index(cx,cy,cz)][0]+
            node_destined_moms[it][im][grid->get_cell_index(cx,cy,iz)][1]+
            node_destined_moms[it][im][grid->get_cell_index(cx,iy,cz)][2]+
            node_destined_moms[it][im][grid->get_cell_index(cx,iy,iz)][3]+
            node_destined_moms[it][im][grid->get_cell_index(ix,cy,cz)][4]+
            node_destined_moms[it][im][grid->get_cell_index(ix,cy,iz)][5]+
            node_destined_moms[it][im][grid->get_cell_index(ix,iy,cz)][6]+
            node_destined_moms[it][im][grid->get_cell_index(ix,iy,iz)][7]);
        }
      }
      else
      {
        // This scatters node data; each threads sums for
        // one moment. To scale beyond 10 threads would need
        // to do even and odd cells separately (to prevent
        // threads from writing to same node).
        //
          for(int it=0;it<num_threads;it++)
          #pragma omp for
          for(int im=0;im<10;im++)
          for(int cx=1;cx<nxc-1;cx++)
          for(int cy=1;cy<nyc-1;cy++)
          for(int cz=1;cz<nzc-1;cz++)
          {
            const int ix=cx+1;
            const int iy=cy+1;
            const int iz=cz+1;
            // scatter data to neighboring nodes
            const int ic = grid->get_cell_index(cx,cy,cz);
            out_moments[im][ix][iy][iz] += invVOL*node_destined_moms[it][im][ic][0];
            out_moments[im][ix][iy][cz] += invVOL*node_destined_moms[it][im][ic][1];
            out_moments[im][ix][cy][iz] += invVOL*node_destined_moms[it][im][ic][2];
            out_moments[im][ix][cy][cz] += invVOL*node_destined_moms[it][im][ic][3];
            out_moments[im][cx][iy][iz] += invVOL*node_destined_moms[it][im][ic][4];
            out_moments[im][cx][iy][cz] += invVOL*node_destined_moms[it][im][ic][5];
            out_moments[im][cx][cy][iz] += invVOL*node_destined_moms[it][im][ic][6];
            out_moments[im][cx][cy][cz] += invVOL*node_destined_moms[it][im][ic][7];
          }
      }
      if(!thread_num) timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
      // uncomment this and remove the loop below
      // when we change to use asynchronous communication.
      // communicateGhostP2G(is, vct);
  }
  #pragma omp master
  delArray4<double>(node_destined_moms);
}
//
// Compare the vectorization notes at the top of mover_PC().
//
// This was Particles3Dcomm::interpP2G()
void SpeciesMoms::sumMoments(const Particles3Dcomm& pcls)
{
  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / grid->getDX();
  const double inv_dy = 1.0 / grid->getDY();
  const double inv_dz = 1.0 / grid->getDZ();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP
  {
    const Particles3Dcomm& pcls = part[i];
    assert_eq(pcls.get_particleType(), ParticleType::SoA);
    const int is = pcls.get_species_num();

    double const*const x = pcls.getXall();
    double const*const y = pcls.getYall();
    double const*const z = pcls.getZall();
    double const*const u = pcls.getUall();
    double const*const v = pcls.getVall();
    double const*const w = pcls.getWall();
    double const*const q = pcls.getQall();

    const int nop = pcls.getNOP();

    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
    int thread_num = omp_get_thread_num();
    Moments10& speciesMoments10 = fetch_moments10Array(thread_num);
    arr4_double moments = speciesMoments10.fetch_arr();
    //
    // moments.setmode(ompmode::mine);
    // moments.setall(0.);
    // 
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    //
    // This barrier is not needed
    #pragma omp barrier
    // The following loop is expensive, so it is wise to assume that the
    // compiler is stupid.  Therefore we should on the one hand
    // expand things out and on the other hand avoid repeating computations.
    #pragma omp for // used nowait with the old way
    for (int i = 0; i < nop; i++)
    {
      // compute the quadratic moments of velocity
      //
      const double ui=u[i];
      const double vi=v[i];
      const double wi=w[i];
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[10];
      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;

      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      const double xi0   = x[i] - grid->getXN(ix-1);
      const double eta0  = y[i] - grid->getYN(iy-1);
      const double zeta0 = z[i] - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - x[i];
      const double eta1  = grid->getYN(iy) - y[i];
      const double zeta1 = grid->getZN(iz) - z[i];
      const double qi = q[i];
      const double invVOLqi = invVOL*qi;
      const double weight0 = invVOLqi * xi0;
      const double weight1 = invVOLqi * xi1;
      const double weight00 = weight0*eta0;
      const double weight01 = weight0*eta1;
      const double weight10 = weight1*eta0;
      const double weight11 = weight1*eta1;
      double weights[8];
      weights[0] = weight00*zeta0; // weight000
      weights[1] = weight00*zeta1; // weight001
      weights[2] = weight01*zeta0; // weight010
      weights[3] = weight01*zeta1; // weight011
      weights[4] = weight10*zeta0; // weight100
      weights[5] = weight10*zeta1; // weight101
      weights[6] = weight11*zeta0; // weight110
      weights[7] = weight11*zeta1; // weight111
      //weights[0] = xi0 * eta0 * zeta0 * qi * invVOL; // weight000
      //weights[1] = xi0 * eta0 * zeta1 * qi * invVOL; // weight001
      //weights[2] = xi0 * eta1 * zeta0 * qi * invVOL; // weight010
      //weights[3] = xi0 * eta1 * zeta1 * qi * invVOL; // weight011
      //weights[4] = xi1 * eta0 * zeta0 * qi * invVOL; // weight100
      //weights[5] = xi1 * eta0 * zeta1 * qi * invVOL; // weight101
      //weights[6] = xi1 * eta1 * zeta0 * qi * invVOL; // weight110
      //weights[7] = xi1 * eta1 * zeta1 * qi * invVOL; // weight111

      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        arr2_double_fetch moments00 = moments[ix  ][iy  ];
        arr2_double_fetch moments01 = moments[ix  ][iy-1];
        arr2_double_fetch moments10 = moments[ix-1][iy  ];
        arr2_double_fetch moments11 = moments[ix-1][iy-1];
        momentsArray[0] = moments00[iz  ]; // moments000 
        momentsArray[1] = moments00[iz-1]; // moments001 
        momentsArray[2] = moments01[iz  ]; // moments010 
        momentsArray[3] = moments01[iz-1]; // moments011 
        momentsArray[4] = moments10[iz  ]; // moments100 
        momentsArray[5] = moments10[iz-1]; // moments101 
        momentsArray[6] = moments11[iz  ]; // moments110 
        momentsArray[7] = moments11[iz-1]; // moments111 

        for(int m=0; m<10; m++)
        for(int c=0; c<8; c++)
        {
          momentsArray[c][m] += velmoments[m]*weights[c];
        }
      }
    }
    #pragma omp master
    timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

    // reduction
    #pragma omp master
    timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

    // reduce moments in parallel
    //
    for(int thread_num=0;thread_num<get_sizeMomentsArray();thread_num++)
    {
      arr4_double moments = fetch_moments10Array(thread_num).fetch_arr();
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++)
      for(int j=0;j<nyn;j++)
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] += invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] += invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9];
      }
    }
    //
    // This was the old way of reducing;
    // did not scale well to large number of threads
    //{
    //  #pragma omp critical (reduceMoment0)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { rhons[is][i][j][k] += invVOL*moments[i][j][k][0]; }}
    //  #pragma omp critical (reduceMoment1)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1]; }}
    //  #pragma omp critical (reduceMoment2)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jys  [is][i][j][k] += invVOL*moments[i][j][k][2]; }}
    //  #pragma omp critical (reduceMoment3)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3]; }}
    //  #pragma omp critical (reduceMoment4)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4]; }}
    //  #pragma omp critical (reduceMoment5)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5]; }}
    //  #pragma omp critical (reduceMoment6)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6]; }}
    //  #pragma omp critical (reduceMoment7)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7]; }}
    //  #pragma omp critical (reduceMoment8)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8]; }}
    //  #pragma omp critical (reduceMoment9)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9]; }}
    //}
    #pragma omp master
    timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
  }
}

void SpeciesMoms::sumMoments_AoS(const Particles3Dcomm& pcls)
{
  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / grid->getDX();
  const double inv_dy = 1.0 / grid->getDY();
  const double inv_dz = 1.0 / grid->getDZ();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  {
    assert_eq(pcls.get_particleType(), ParticleType::AoS);
    const int is = pcls.get_species_num();

    const int nop = pcls.getNOP();

    int thread_num = omp_get_thread_num();
    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
    Moments10& speciesMoments10 = fetch_moments10Array(thread_num);
    arr4_double moments = speciesMoments10.fetch_arr();
    //
    // moments.setmode(ompmode::mine);
    // moments.setall(0.);
    // 
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    //
    #pragma omp barrier
    #pragma omp for
    for (int pidx = 0; pidx < nop; pidx++)
    {
      const SpeciesParticle& pcl = pcls.get_pcl(pidx);
      // compute the quadratic moments of velocity
      //
      const double ui=pcl.get_u();
      const double vi=pcl.get_v();
      const double wi=pcl.get_w();
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[10];
      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;

      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((pcl.get_x() - xstart) * inv_dx));
      const int iy = 2 + int (floor((pcl.get_y() - ystart) * inv_dy));
      const int iz = 2 + int (floor((pcl.get_z() - zstart) * inv_dz));
      const double xi0   = pcl.get_x() - grid->getXN(ix-1);
      const double eta0  = pcl.get_y() - grid->getYN(iy-1);
      const double zeta0 = pcl.get_z() - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - pcl.get_x();
      const double eta1  = grid->getYN(iy) - pcl.get_y();
      const double zeta1 = grid->getZN(iz) - pcl.get_z();
      const double qi = pcl.get_q();
      const double invVOLqi = invVOL*qi;
      const double weight0 = invVOLqi * xi0;
      const double weight1 = invVOLqi * xi1;
      const double weight00 = weight0*eta0;
      const double weight01 = weight0*eta1;
      const double weight10 = weight1*eta0;
      const double weight11 = weight1*eta1;
      double weights[8];
      weights[0] = weight00*zeta0; // weight000
      weights[1] = weight00*zeta1; // weight001
      weights[2] = weight01*zeta0; // weight010
      weights[3] = weight01*zeta1; // weight011
      weights[4] = weight10*zeta0; // weight100
      weights[5] = weight10*zeta1; // weight101
      weights[6] = weight11*zeta0; // weight110
      weights[7] = weight11*zeta1; // weight111
      //weights[0] = xi0 * eta0 * zeta0 * qi * invVOL; // weight000
      //weights[1] = xi0 * eta0 * zeta1 * qi * invVOL; // weight001
      //weights[2] = xi0 * eta1 * zeta0 * qi * invVOL; // weight010
      //weights[3] = xi0 * eta1 * zeta1 * qi * invVOL; // weight011
      //weights[4] = xi1 * eta0 * zeta0 * qi * invVOL; // weight100
      //weights[5] = xi1 * eta0 * zeta1 * qi * invVOL; // weight101
      //weights[6] = xi1 * eta1 * zeta0 * qi * invVOL; // weight110
      //weights[7] = xi1 * eta1 * zeta1 * qi * invVOL; // weight111

      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        arr2_double_fetch moments00 = moments[ix  ][iy  ];
        arr2_double_fetch moments01 = moments[ix  ][iy-1];
        arr2_double_fetch moments10 = moments[ix-1][iy  ];
        arr2_double_fetch moments11 = moments[ix-1][iy-1];
        momentsArray[0] = moments00[iz  ]; // moments000 
        momentsArray[1] = moments00[iz-1]; // moments001 
        momentsArray[2] = moments01[iz  ]; // moments010 
        momentsArray[3] = moments01[iz-1]; // moments011 
        momentsArray[4] = moments10[iz  ]; // moments100 
        momentsArray[5] = moments10[iz-1]; // moments101 
        momentsArray[6] = moments11[iz  ]; // moments110 
        momentsArray[7] = moments11[iz-1]; // moments111 

        for(int m=0; m<10; m++)
        for(int c=0; c<8; c++)
        {
          momentsArray[c][m] += velmoments[m]*weights[c];
        }
      }
    }
    if(!thread_num) timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

    // reduction
    if(!thread_num) timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

    // reduce moments in parallel
    //
    for(int thread_num=0;thread_num<get_sizeMomentsArray();thread_num++)
    {
      arr4_double moments = fetch_moments10Array(thread_num).fetch_arr();
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++)
      for(int j=0;j<nyn;j++)
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] += invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] += invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9];
      }
    }
    if(!thread_num) timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
  }
}

#ifdef __MIC__
// add moment weights to all ten moments for the cell of the particle
// (assumes that particle data is aligned with cache boundary and
// begins with the velocity components)
inline void addto_cell_moments(
  F64vec8* cell_moments,
  F64vec8 weights,
  F64vec8 vel)
{
  // broadcast particle velocities
  const F64vec8 u = F64vec8(vel[0]);
  const F64vec8 v = F64vec8(vel[1]);
  const F64vec8 w = F64vec8(vel[2]);
  // construct kronecker product of moments and weights
  const F64vec8 u_weights = u*weights;
  const F64vec8 v_weights = v*weights;
  const F64vec8 w_weights = w*weights;
  const F64vec8 uu_weights = u*u_weights;
  const F64vec8 uv_weights = u*v_weights;
  const F64vec8 uw_weights = u*w_weights;
  const F64vec8 vv_weights = v*v_weights;
  const F64vec8 vw_weights = v*w_weights;
  const F64vec8 ww_weights = w*w_weights;
  // add moment weights to accumulated moment weights in mesh mesh
  cell_moments[0] += weights;
  cell_moments[1] += u_weights;
  cell_moments[2] += v_weights;
  cell_moments[3] += w_weights;
  cell_moments[4] += uu_weights;
  cell_moments[5] += uv_weights;
  cell_moments[6] += uw_weights;
  cell_moments[7] += vv_weights;
  cell_moments[8] += vw_weights;
  cell_moments[9] += ww_weights;
}
#endif // __MIC__

// sum moments of AoS using MIC intrinsics
// 
// We could rewrite this without intrinsics also.  The core idea
// of this algorithm is that instead of scattering the data of
// each particle to its nodes, in each cell we accumulate the
// data that would be scattered and then scatter it at the end.
// By waiting to scatter, with each particle we work with an
// aligned 10x8 matrix rather than a 8x10 matrix, which means
// that for each particle we make 10 vector stores rather than
// 8*2=16 or 8*3=24 vector stores (for unaligned data).  This
// also avoids the expense of computing node indices for each
// particle.
//
// 1. compute vector of 8 weights using position
// 2. form kronecker product of weights with moments
//    by scaling the weights by each velocity moment;
//    add each to accumulated weights for this cell
// 3. after sum is complete, transpose weight-moment
//    product in each cell and distribute to its 8 nodes.
//    An optimized way:
//    A. transpose the first 8 weighted moments with fast 8x8
//       matrix transpose.
//    B. transpose 2x8 matrix of the last two weighted moments
//       and then use 8 masked vector adds to accumulate
//       to weights at nodes.
//    But the optimized way might be overkill since distributing
//    the sums from the cells to the nodes should not dominate
//    if the number of particles per mesh cell is large;
//    if the number of particles per mesh cell is small,
//    then a fully vectorized moment sum is hard to justify anyway.
//
// See notes at the top of sumMoments().
//
void SpeciesMoms::sumMoments_AoS_intr(const Particles3Dcomm& pcls)
{
#ifndef __MIC__
  eprintf("not implemented");
#else
  const Grid *grid = &get_grid();

  // define global parameters
  //
  const double inv_dx = 1.0 / grid->getDX();
  const double inv_dy = 1.0 / grid->getDY();
  const double inv_dz = 1.0 / grid->getDZ();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // Here and below x stands for all 3 physical position coordinates
  const F64vec8 dx_inv = make_F64vec8(inv_dx, inv_dy, inv_dz);
  // starting physical position of proper subdomain ("pdom", without ghosts)
  const F64vec8 pdom_xlow = make_F64vec8(xstart,ystart, zstart);
  //
  // X = canonical coordinates.
  //
  // starting position of cell in lower corner
  // of proper subdomain (without ghosts);
  // probably this is an integer value, but we won't rely on it.
  const F64vec8 pdom_Xlow = dx_inv*pdom_xlow;
  // g = including ghosts
  // starting position of cell in low corner
  const F64vec8 gdom_Xlow = pdom_Xlow - F64vec8(1.);
  // starting position of cell in high corner of physical domain
  // in canonical coordinates of ghost domain
  const F64vec8 nXcm1 = make_F64vec8(nxc-1,nyc-1,nzc-1);

  // allocate memory per mesh cell for accumulating moments
  //
  const int num_threads = omp_get_max_threads();
  array4<F64vec8>* cell_moments_per_thr;
  #pragma omp single
  {
    cell_moments_per_thr
    = (array4<F64vec8>*) malloc(num_threads*sizeof(array4<F64vec8>));
  }
  #pragma omp single //#pragma omp for // (is memory allocation thread-safe?)
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // use placement new to allocate array to accumulate moments for thread
    new(&cell_moments_per_thr[thread_num]) array4<F64vec8>(nxc,nyc,nzc,10);
  }
  //
  // allocate memory per mesh node for accumulating moments
  //
  array3<F64vec8>* node_moments_first8_per_thr;
  array4<double>* node_moments_last2_per_thr;
  #pragma omp single
  {
    node_moments_first8_per_thr
    = (array3<F64vec8>*) malloc(num_threads*sizeof(array3<F64vec8>));
    node_moments_last2_per_thr
    = (array4<double>*) malloc(num_threads*sizeof(array4<double>));
  }
  #pragma omp single //#pragma omp for // (is memory allocation thread-safe?)
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // use placement new to allocate array to accumulate moments for thread
    new(&node_moments_first8_per_thr[thread_num]) array3<F64vec8>(nxn,nyn,nzn);
    new(&node_moments_last2_per_thr[thread_num]) array4<double>(nxn,nyn,nzn,2);
  }

  // The moments of a particle must be distributed to the 8 nodes of the cell
  // in proportion to the weight of each node.
  //
  // Refer to the kronecker product of weights and moments as
  // "weighted moments" or "moment weights".
  //
  // Each thread accumulates moment weights in cells.
  //
  // Because particles are not assumed to be sorted by mesh cell,
  // we have to wait until all particles have been processed
  // before we transpose moment weights to weighted moments;
  // the memory that we must allocate to sum moments is thus
  // num_thread*8 times as much as if particles were pre-sorted
  // by mesh cell (and num_threads times as much as if particles
  // were sorted by thread subdomain).
  //
  // #assert omp parallel
  {
    // array4<F64vec8> cell_moments(nxc,nyc,nzc,10);
    const int this_thread = omp_get_thread_num();
    assert_lt(this_thread,num_threads);
    array4<F64vec8>& cell_moments = cell_moments_per_thr[this_thread];

    {
      assert_eq(pcls.get_particleType(), ParticleType::AoS);
      const int is = pcls.get_species_num();

      // moments.setmode(ompmode::mine);
      // moments.setall(0.);
      // 
      F64vec8 *cell_moments1d = &cell_moments[0][0][0][0];
      int moments1dsize = cell_moments.get_size();
      for(int i=0; i<moments1dsize; i++) cell_moments1d[i]=F64vec8(0.);
      //
      // number or particles processed at a time
      const int num_pcls_per_loop = 2;
      const vector_SpeciesParticle& pcl_list = pcls.get_pcl_list();
      const int nop = pcl_list.size();
      // if the number of particles is odd, then make
      // sure that the data after the last particle
      // will not contribute to the moments.
      #pragma omp single // the implied omp barrier is needed
      {
        // make sure that we will not overrun the array
        assert_divides(num_pcls_per_loop,pcl_list.capacity());
        // round up number of particles
        int nop_rounded_up = roundup_to_multiple(nop,num_pcls_per_loop);
        for(int pidx=nop; pidx<nop_rounded_up; pidx++)
        {
          // (This is a benign violation of particle
          // encapsulation and requires a cast).
          SpeciesParticle& pcl = (SpeciesParticle&) pcl_list[pidx];
          pcl.set_to_zero();
        }
      }
      #pragma omp for
      for (int pidx = 0; pidx < nop; pidx+=2)
      {
        // cast particles as vectors
        // (assumes each particle exactly fits a cache line)
        const F64vec8& pcl0 = (const F64vec8&)pcl_list[pidx];
        const F64vec8& pcl1 = (const F64vec8&)pcl_list[pidx+1];
        // gather position data from particles
        // (assumes position vectors are in upper half)
        const F64vec8 xpos = cat_hgh_halves(pcl0,pcl1);

        // convert to canonical coordinates relative to subdomain with ghosts
        const F64vec8 gX = dx_inv*xpos - gdom_Xlow;
        F64vec8 cellXstart = floor(gX);
        // all particles at this point should be inside the
        // proper subdomain of this process, but maybe we
        // will need to enforce this because of inconsistency
        // of floating point arithmetic?
        //cellXstart = maximum(cellXstart,F64vec8(1.));
        //cellXstart = minimum(cellXstart,nXcm1);
        assert(!test_lt(cellXstart,F64vec8(1.)));
        assert(!test_gt(cellXstart,nXcm1));

        // get weights for field_components based on particle position
        //
        F64vec8 weights[2];
        const F64vec8 X = gX - cellXstart;
        construct_weights_for_2pcls(weights, X);

        // add scaled weights to all ten moments for the cell of each particle
        //
        // the cell that we will write to
        const I32vec16 cell = round_to_nearest(cellXstart);
        const int* c=(int*)&cell;
        F64vec8* cell_moments0 = &cell_moments[c[0]][c[1]][c[2]][0];
        F64vec8* cell_moments1 = &cell_moments[c[4]][c[5]][c[6]][0];
        addto_cell_moments(cell_moments0, weights[0], pcl0);
        addto_cell_moments(cell_moments1, weights[1], pcl1);
      }
      if(!this_thread) timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

      // reduction
      if(!this_thread) timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

      // reduce moments in parallel
      //
      // this code currently makes no sense for multiple threads.
      assert_eq(num_threads,1);
      {
        // For each thread, distribute moments from cells to nodes
        // and then sum moments at each node over all threads.
        //
        // (Alternatively we could sum over all threads and then
        // distribute to nodes; this alternative would be preferable
        // for vectorization efficiency but more difficult to parallelize
        // across threads).

        // initialize moment accumulators
        //
        memset(&node_moments_first8_per_thr[this_thread][0][0][0],
          0, sizeof(F64vec8)*node_moments_first8_per_thr[0].get_size());
        memset(&node_moments_last2_per_thr[this_thread][0][0][0][0],
          0, sizeof(double)*node_moments_last2_per_thr[0].get_size());

        // distribute moments from cells to nodes
        //
        #pragma omp for collapse(2)
        for(int cx=1;cx<nxc;cx++)
        for(int cy=1;cy<nyc;cy++)
        for(int cz=1;cz<nzc;cz++)
        {
          const int ix=cx+1;
          const int iy=cy+1;
          const int iz=cz+1;
          F64vec8* cell_mom = &cell_moments[cx][cy][cz][0];

          // scatter the cell's first 8 moments to its nodes
          // for each thread
          {
            F64vec8* cell_mom_first8 = cell_mom;
            // regard cell_mom_first8 as a pointer to 8x8 data and transpose
            transpose_8x8_double((double(*)[8]) cell_mom_first8);
            // scatter the moment vectors to the nodes
            array3<F64vec8>& node_moments_first8 = node_moments_first8_per_thr[this_thread];
            arr_fetch2(F64vec8) node_moments0 = node_moments_first8[ix];
            arr_fetch2(F64vec8) node_moments1 = node_moments_first8[cx];
            arr_fetch1(F64vec8) node_moments00 = node_moments0[iy];
            arr_fetch1(F64vec8) node_moments01 = node_moments0[cy];
            arr_fetch1(F64vec8) node_moments10 = node_moments1[iy];
            arr_fetch1(F64vec8) node_moments11 = node_moments1[cy];
            node_moments00[iz] += cell_mom_first8[0]; // node_moments_first8[ix][iy][iz]
            node_moments00[cz] += cell_mom_first8[1]; // node_moments_first8[ix][iy][cz]
            node_moments01[iz] += cell_mom_first8[2]; // node_moments_first8[ix][cy][iz]
            node_moments01[cz] += cell_mom_first8[3]; // node_moments_first8[ix][cy][cz]
            node_moments10[iz] += cell_mom_first8[4]; // node_moments_first8[cx][iy][iz]
            node_moments10[cz] += cell_mom_first8[5]; // node_moments_first8[cx][iy][cz]
            node_moments11[iz] += cell_mom_first8[6]; // node_moments_first8[cx][cy][iz]
            node_moments11[cz] += cell_mom_first8[7]; // node_moments_first8[cx][cy][cz]
          }

          // scatter the cell's last 2 moments to its nodes
          {
            array4<double>& node_moments_last2 = node_moments_last2_per_thr[this_thread];
            arr3_double_fetch node_moments0 = node_moments_last2[ix];
            arr3_double_fetch node_moments1 = node_moments_last2[cx];
            arr2_double_fetch node_moments00 = node_moments0[iy];
            arr2_double_fetch node_moments01 = node_moments0[cy];
            arr2_double_fetch node_moments10 = node_moments1[iy];
            arr2_double_fetch node_moments11 = node_moments1[cy];
            double* node_moments000 = node_moments00[iz];
            double* node_moments001 = node_moments00[cz];
            double* node_moments010 = node_moments01[iz];
            double* node_moments011 = node_moments01[cz];
            double* node_moments100 = node_moments10[iz];
            double* node_moments101 = node_moments10[cz];
            double* node_moments110 = node_moments11[iz];
            double* node_moments111 = node_moments11[cz];

            const F64vec8 mom8 = cell_mom[8];
            const F64vec8 mom9 = cell_mom[9];

            bool naive_last2 = true;
            if(naive_last2)
            {
              node_moments000[0] += mom8[0]; node_moments000[1] += mom9[0];
              node_moments001[0] += mom8[1]; node_moments001[1] += mom9[1];
              node_moments010[0] += mom8[2]; node_moments010[1] += mom9[2];
              node_moments011[0] += mom8[3]; node_moments011[1] += mom9[3];
              node_moments100[0] += mom8[4]; node_moments100[1] += mom9[4];
              node_moments101[0] += mom8[5]; node_moments101[1] += mom9[5];
              node_moments110[0] += mom8[6]; node_moments110[1] += mom9[6];
              node_moments111[0] += mom8[7]; node_moments111[1] += mom9[7];
            }
            else
            {
              // Let a=moment#8 and b=moment#9.
              // Number the nodes 0 through 7.
              //
              // This transpose changes data from the form
              //   [a0 a1 a2 a3 a4 a5 a6 a7]=mom8
              //   [b0 b1 b2 b3 b4 b5 b6 b7]=mom9
              // into the form
              //   [a0 b0 a2 b2 a4 b4 a6 b6]=out8
              //   [a1 b1 a3 b3 a5 b5 a7 b7]=out9
              F64vec8 out8, out9;
              trans2x2(out8, out9, mom8, mom9);

              // probably the compiler is not smart enough to recognize that
              // each line can be done with a single vector instruction:
              node_moments000[0] += out8[0]; node_moments000[1] += out8[1];
              node_moments001[0] += out9[0]; node_moments001[1] += out9[1];
              node_moments010[0] += out8[2]; node_moments010[1] += out8[3];
              node_moments011[0] += out9[2]; node_moments011[1] += out9[3];
              node_moments100[0] += out8[4]; node_moments100[1] += out8[5];
              node_moments101[0] += out9[4]; node_moments101[1] += out9[5];
              node_moments110[0] += out8[6]; node_moments110[1] += out8[7];
              node_moments111[0] += out9[6]; node_moments111[1] += out9[7];
            }
          }
        }

        // at each node add moments to moments of first thread
        //
        #pragma omp for collapse(2)
        for(int nx=1;nx<nxn;nx++)
        for(int ny=1;ny<nyn;ny++)
        {
          arr_fetch1(F64vec8) node_moments8_for_master
            = node_moments_first8_per_thr[0][nx][ny];
          arr_fetch2(double) node_moments2_for_master
            = node_moments_last2_per_thr[0][nx][ny];
          for(int thread_num=1;thread_num<num_threads;thread_num++)
          {
            arr_fetch1(F64vec8) node_moments8_for_thr
              = node_moments_first8_per_thr[thread_num][nx][ny];
            arr_fetch2(double) node_moments2_for_thr
              = node_moments_last2_per_thr[thread_num][nx][ny];
            for(int nz=1;nz<nzn;nz++)
            {
              node_moments8_for_master[nz] += node_moments8_for_thr[nz];
              node_moments2_for_master[nz][0] += node_moments2_for_thr[nz][0];
              node_moments2_for_master[nz][1] += node_moments2_for_thr[nz][1];
            }
          }
        }

        // transpose moments for field solver
        //
        #pragma omp for collapse(2)
        for(int nx=1;nx<nxn;nx++)
        for(int ny=1;ny<nyn;ny++)
        {
          arr_fetch1(F64vec8) node_moments8_for_master
            = node_moments_first8_per_thr[0][nx][ny];
          arr_fetch2(double) node_moments2_for_master
            = node_moments_last2_per_thr[0][nx][ny];
          arr_fetch1(double) rho_sxy = rhons[is][nx][ny];
          arr_fetch1(double) Jx__sxy = Jxs  [is][nx][ny];
          arr_fetch1(double) Jy__sxy = Jys  [is][nx][ny];
          arr_fetch1(double) Jz__sxy = Jzs  [is][nx][ny];
          arr_fetch1(double) pXX_sxy = pXXsn[is][nx][ny];
          arr_fetch1(double) pXY_sxy = pXYsn[is][nx][ny];
          arr_fetch1(double) pXZ_sxy = pXZsn[is][nx][ny];
          arr_fetch1(double) pYY_sxy = pYYsn[is][nx][ny];
          arr_fetch1(double) pYZ_sxy = pYZsn[is][nx][ny];
          arr_fetch1(double) pZZ_sxy = pZZsn[is][nx][ny];
          for(int nz=0;nz<nzn;nz++)
          {
            rho_sxy[nz] = invVOL*node_moments8_for_master[nz][0];
            Jx__sxy[nz] = invVOL*node_moments8_for_master[nz][1];
            Jy__sxy[nz] = invVOL*node_moments8_for_master[nz][2];
            Jz__sxy[nz] = invVOL*node_moments8_for_master[nz][3];
            pXX_sxy[nz] = invVOL*node_moments8_for_master[nz][4];
            pXY_sxy[nz] = invVOL*node_moments8_for_master[nz][5];
            pXZ_sxy[nz] = invVOL*node_moments8_for_master[nz][6];
            pYY_sxy[nz] = invVOL*node_moments8_for_master[nz][7];
            pYZ_sxy[nz] = invVOL*node_moments2_for_master[nz][0];
            pZZ_sxy[nz] = invVOL*node_moments2_for_master[nz][1];
          }
        }
      }
      if(!this_thread) timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
    }
  }

  #pragma omp single
  {
    // deallocate memory per mesh node for accumulating moments
    //
    for(int thread_num=0;thread_num<num_threads;thread_num++)
    {
      // call destructor to deallocate arrays
      node_moments_first8_per_thr[thread_num].~array3<F64vec8>();
      node_moments_last2_per_thr[thread_num].~array4<double>();
    }
    free(node_moments_first8_per_thr);
    free(node_moments_last2_per_thr);

    // deallocate memory for accumulating moments
    //
    for(int thread_num=0;thread_num<num_threads;thread_num++)
    {
      // deallocate array to accumulate moments for thread
      cell_moments_per_thr[thread_num].~array4<F64vec8>();
    }
    free(cell_moments_per_thr);
  }
#endif // __MIC__
}

inline void compute_moments(double velmoments[10], double weights[8],
  int i,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  ALIGNED(x);
  ALIGNED(y);
  ALIGNED(z);
  ALIGNED(u);
  ALIGNED(v);
  ALIGNED(w);
  ALIGNED(q);
  // compute the quadratic moments of velocity
  //
  const double ui=u[i];
  const double vi=v[i];
  const double wi=w[i];
  const double uui=ui*ui;
  const double uvi=ui*vi;
  const double uwi=ui*wi;
  const double vvi=vi*vi;
  const double vwi=vi*wi;
  const double wwi=wi*wi;
  //double velmoments[10];
  velmoments[0] = 1.;
  velmoments[1] = ui;
  velmoments[2] = vi;
  velmoments[3] = wi;
  velmoments[4] = uui;
  velmoments[5] = uvi;
  velmoments[6] = uwi;
  velmoments[7] = vvi;
  velmoments[8] = vwi;
  velmoments[9] = wwi;

  // compute the weights to distribute the moments
  //
  //double weights[8];
  const double abs_xpos = x[i];
  const double abs_ypos = y[i];
  const double abs_zpos = z[i];
  const double rel_xpos = abs_xpos - xstart;
  const double rel_ypos = abs_ypos - ystart;
  const double rel_zpos = abs_zpos - zstart;
  const double cxm1_pos = rel_xpos * inv_dx;
  const double cym1_pos = rel_ypos * inv_dy;
  const double czm1_pos = rel_zpos * inv_dz;
  //if(true)
  //{
  //  const int cx_inf = int(floor(cxm1_pos));
  //  const int cy_inf = int(floor(cym1_pos));
  //  const int cz_inf = int(floor(czm1_pos));
  //  assert_eq(cx-1,cx_inf);
  //  assert_eq(cy-1,cy_inf);
  //  assert_eq(cz-1,cz_inf);
  //}
  // fraction of the distance from the right of the cell
  const double w1x = cx - cxm1_pos;
  const double w1y = cy - cym1_pos;
  const double w1z = cz - czm1_pos;
  // fraction of distance from the left
  const double w0x = 1-w1x;
  const double w0y = 1-w1y;
  const double w0z = 1-w1z;
  // we are calculating a charge moment.
  const double qi=q[i];
  const double weight0 = qi*w0x;
  const double weight1 = qi*w1x;
  const double weight00 = weight0*w0y;
  const double weight01 = weight0*w1y;
  const double weight10 = weight1*w0y;
  const double weight11 = weight1*w1y;
  weights[0] = weight00*w0z; // weight000
  weights[1] = weight00*w1z; // weight001
  weights[2] = weight01*w0z; // weight010
  weights[3] = weight01*w1z; // weight011
  weights[4] = weight10*w0z; // weight100
  weights[5] = weight10*w1z; // weight101
  weights[6] = weight11*w0z; // weight110
  weights[7] = weight11*w1z; // weight111
}

// add particle to moments
inline void add_moments_for_pcl(double momentsAcc[8][10],
  int i,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  double velmoments[10];
  double weights[8];
  compute_moments(velmoments,weights,
    i, x, y, z, u, v, w, q,
    xstart, ystart, zstart,
    inv_dx, inv_dy, inv_dz,
    cx, cy, cz);

  // add moments for this particle
  {
    // which is the superior order for the following loop?
    for(int c=0; c<8; c++)
    for(int m=0; m<10; m++)
    {
      momentsAcc[c][m] += velmoments[m]*weights[c];
    }
  }
}


// vectorized version of previous method
// 
inline void add_moments_for_pcl_vec(double momentsAccVec[8][10][8],
  double velmoments[10][8], double weights[8][8],
  int i,
  int imod,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  ALIGNED(x);
  ALIGNED(y);
  ALIGNED(z);
  ALIGNED(u);
  ALIGNED(v);
  ALIGNED(w);
  ALIGNED(q);
  // compute the quadratic moments of velocity
  //
  const double ui=u[i];
  const double vi=v[i];
  const double wi=w[i];
  const double uui=ui*ui;
  const double uvi=ui*vi;
  const double uwi=ui*wi;
  const double vvi=vi*vi;
  const double vwi=vi*wi;
  const double wwi=wi*wi;
  //double velmoments[10];
  velmoments[0][imod] = 1.;
  velmoments[1][imod] = ui;
  velmoments[2][imod] = vi;
  velmoments[3][imod] = wi;
  velmoments[4][imod] = uui;
  velmoments[5][imod] = uvi;
  velmoments[6][imod] = uwi;
  velmoments[7][imod] = vvi;
  velmoments[8][imod] = vwi;
  velmoments[9][imod] = wwi;

  // compute the weights to distribute the moments
  //
  //double weights[8];
  const double abs_xpos = x[i];
  const double abs_ypos = y[i];
  const double abs_zpos = z[i];
  const double rel_xpos = abs_xpos - xstart;
  const double rel_ypos = abs_ypos - ystart;
  const double rel_zpos = abs_zpos - zstart;
  const double cxm1_pos = rel_xpos * inv_dx;
  const double cym1_pos = rel_ypos * inv_dy;
  const double czm1_pos = rel_zpos * inv_dz;
  //if(true)
  //{
  //  const int cx_inf = int(floor(cxm1_pos));
  //  const int cy_inf = int(floor(cym1_pos));
  //  const int cz_inf = int(floor(czm1_pos));
  //  assert_eq(cx-1,cx_inf);
  //  assert_eq(cy-1,cy_inf);
  //  assert_eq(cz-1,cz_inf);
  //}
  // fraction of the distance from the right of the cell
  const double w1x = cx - cxm1_pos;
  const double w1y = cy - cym1_pos;
  const double w1z = cz - czm1_pos;
  // fraction of distance from the left
  const double w0x = 1-w1x;
  const double w0y = 1-w1y;
  const double w0z = 1-w1z;
  // we are calculating a charge moment.
  const double qi=q[i];
  const double weight0 = qi*w0x;
  const double weight1 = qi*w1x;
  const double weight00 = weight0*w0y;
  const double weight01 = weight0*w1y;
  const double weight10 = weight1*w0y;
  const double weight11 = weight1*w1y;
  weights[0][imod] = weight00*w0z; // weight000
  weights[1][imod] = weight00*w1z; // weight001
  weights[2][imod] = weight01*w0z; // weight010
  weights[3][imod] = weight01*w1z; // weight011
  weights[4][imod] = weight10*w0z; // weight100
  weights[5][imod] = weight10*w1z; // weight101
  weights[6][imod] = weight11*w0z; // weight110
  weights[7][imod] = weight11*w1z; // weight111

  // add moments for this particle
  {
    for(int c=0; c<8; c++)
    for(int m=0; m<10; m++)
    {
      momentsAccVec[c][m][imod] += velmoments[m][imod]*weights[c][imod];
    }
  }
}

//void SpeciesMoms::sumMoments_vectorized(const Particles3Dcomm* part)
//{
//  const Grid *grid = &get_grid();
//
//  const double inv_dx = grid->get_invdx();
//  const double inv_dy = grid->get_invdy();
//  const double inv_dz = grid->get_invdz();
//  const int nxn = grid->getNXN();
//  const int nyn = grid->getNYN();
//  const int nzn = grid->getNZN();
//  const double xstart = grid->getXstart();
//  const double ystart = grid->getYstart();
//  const double zstart = grid->getZstart();
//  #pragma omp parallel
//  {
//  for (int species_idx = 0; species_idx < ns; species_idx++)
//  {
//    const Particles3Dcomm& pcls = part[species_idx];
//    assert_eq(pcls.get_particleType(), ParticleType::SoA);
//    const int is = pcls.get_species_num();
//    assert_eq(species_idx,is);
//
//    double const*const x = pcls.getXall();
//    double const*const y = pcls.getYall();
//    double const*const z = pcls.getZall();
//    double const*const u = pcls.getUall();
//    double const*const v = pcls.getVall();
//    double const*const w = pcls.getWall();
//    double const*const q = pcls.getQall();
//
//    const int nop = pcls.getNOP();
//    #pragma omp master
//    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
//    Moments10& speciesMoments10 = fetch_moments10Array(0);
//    arr4_double moments = speciesMoments10.fetch_arr();
//    //
//    // moments.setmode(ompmode::ompfor);
//    //moments.setall(0.);
//    double *moments1d = &moments[0][0][0][0];
//    int moments1dsize = moments.get_size();
//    #pragma omp for // because shared
//    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
//    
//    // prevent threads from writing to the same location
//    for(int cxmod2=0; cxmod2<2; cxmod2++)
//    for(int cymod2=0; cymod2<2; cymod2++)
//    // each mesh cell is handled by its own thread
//    #pragma omp for collapse(2)
//    for(int cx=cxmod2;cx<nxc;cx+=2)
//    for(int cy=cymod2;cy<nyc;cy+=2)
//    for(int cz=0;cz<nzc;cz++)
//    {
//     //dprint(cz);
//     // index of interface to right of cell
//     const int ix = cx + 1;
//     const int iy = cy + 1;
//     const int iz = cz + 1;
//     {
//      // reference the 8 nodes to which we will
//      // write moment data for particles in this mesh cell.
//      //
//      arr1_double_fetch momentsArray[8];
//      arr2_double_fetch moments00 = moments[ix][iy];
//      arr2_double_fetch moments01 = moments[ix][cy];
//      arr2_double_fetch moments10 = moments[cx][iy];
//      arr2_double_fetch moments11 = moments[cx][cy];
//      momentsArray[0] = moments00[iz]; // moments000 
//      momentsArray[1] = moments00[cz]; // moments001 
//      momentsArray[2] = moments01[iz]; // moments010 
//      momentsArray[3] = moments01[cz]; // moments011 
//      momentsArray[4] = moments10[iz]; // moments100 
//      momentsArray[5] = moments10[cz]; // moments101 
//      momentsArray[6] = moments11[iz]; // moments110 
//      momentsArray[7] = moments11[cz]; // moments111 
//
//      const int numpcls_in_cell = pcls.get_numpcls_in_bucket(cx,cy,cz);
//      const int bucket_offset = pcls.get_bucket_offset(cx,cy,cz);
//      const int bucket_end = bucket_offset+numpcls_in_cell;
//
//      bool vectorized=false;
//      if(!vectorized)
//      {
//        // accumulators for moments per each of 8 threads
//        double momentsAcc[8][10];
//        memset(momentsAcc,0,sizeof(double)*8*10);
//        for(int i=bucket_offset; i<bucket_end; i++)
//        {
//          add_moments_for_pcl(momentsAcc, i,
//            x, y, z, u, v, w, q,
//            xstart, ystart, zstart,
//            inv_dx, inv_dy, inv_dz,
//            cx, cy, cz);
//        }
//        for(int c=0; c<8; c++)
//        for(int m=0; m<10; m++)
//        {
//          momentsArray[c][m] += momentsAcc[c][m];
//        }
//      }
//      if(vectorized)
//      {
//        double velmoments[10][8];
//        double weights[8][8];
//        double momentsAccVec[8][10][8];
//        memset(momentsAccVec,0,sizeof(double)*8*10*8);
//        #pragma simd
//        for(int i=bucket_offset; i<bucket_end; i++)
//        {
//          add_moments_for_pcl_vec(momentsAccVec, velmoments, weights,
//            i, i%8,
//            x, y, z, u, v, w, q,
//            xstart, ystart, zstart,
//            inv_dx, inv_dy, inv_dz,
//            cx, cy, cz);
//        }
//        for(int c=0; c<8; c++)
//        for(int m=0; m<10; m++)
//        for(int i=0; i<8; i++)
//        {
//          momentsArray[c][m] += momentsAccVec[c][m][i];
//        }
//      }
//     }
//    }
//    #pragma omp master
//    { timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION); }
//
//    // reduction
//    #pragma omp master
//    { timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION); }
//    {
//      #pragma omp for collapse(2)
//      for(int i=0;i<nxn;i++){
//      for(int j=0;j<nyn;j++){
//      for(int k=0;k<nzn;k++)
//      {
//        rhons[is][i][j][k] = invVOL*moments[i][j][k][0];
//        Jxs  [is][i][j][k] = invVOL*moments[i][j][k][1];
//        Jys  [is][i][j][k] = invVOL*moments[i][j][k][2];
//        Jzs  [is][i][j][k] = invVOL*moments[i][j][k][3];
//        pXXsn[is][i][j][k] = invVOL*moments[i][j][k][4];
//        pXYsn[is][i][j][k] = invVOL*moments[i][j][k][5];
//        pXZsn[is][i][j][k] = invVOL*moments[i][j][k][6];
//        pYYsn[is][i][j][k] = invVOL*moments[i][j][k][7];
//        pYZsn[is][i][j][k] = invVOL*moments[i][j][k][8];
//        pZZsn[is][i][j][k] = invVOL*moments[i][j][k][9];
//      }}}
//    }
//    #pragma omp master
//    { timeTasks_end_task(TimeTasks::MOMENT_REDUCTION); }
//    // uncomment this and remove the loop below
//    // when we change to use asynchronous communication.
//    // communicateGhostP2G(is);
//  }
//  }
//  for (int i = 0; i < ns; i++)
//  {
//    communicateGhostP2G(i);
//  }
//}
//
//void SpeciesMoms::sumMoments_vectorized_AoS(const Particles3Dcomm* part)
//{
//  const Grid *grid = &get_grid();
//
//  const double inv_dx = grid->get_invdx();
//  const double inv_dy = grid->get_invdy();
//  const double inv_dz = grid->get_invdz();
//  const int nxn = grid->getNXN();
//  const int nyn = grid->getNYN();
//  const int nzn = grid->getNZN();
//  const double xstart = grid->getXstart();
//  const double ystart = grid->getYstart();
//  const double zstart = grid->getZstart();
//  #pragma omp parallel
//  {
//  for (int species_idx = 0; species_idx < ns; species_idx++)
//  {
//    const Particles3Dcomm& pcls = part[species_idx];
//    assert_eq(pcls.get_particleType(), ParticleType::AoS);
//    const int is = pcls.get_species_num();
//    assert_eq(species_idx,is);
//
//    const int nop = pcls.getNOP();
//    #pragma omp master
//    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
//    Moments10& speciesMoments10 = fetch_moments10Array(0);
//    arr4_double moments = speciesMoments10.fetch_arr();
//    //
//    // moments.setmode(ompmode::ompfor);
//    //moments.setall(0.);
//    double *moments1d = &moments[0][0][0][0];
//    int moments1dsize = moments.get_size();
//    #pragma omp for // because shared
//    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
//    
//    // prevent threads from writing to the same location
//    for(int cxmod2=0; cxmod2<2; cxmod2++)
//    for(int cymod2=0; cymod2<2; cymod2++)
//    // each mesh cell is handled by its own thread
//    #pragma omp for collapse(2)
//    for(int cx=cxmod2;cx<nxc;cx+=2)
//    for(int cy=cymod2;cy<nyc;cy+=2)
//    for(int cz=0;cz<nzc;cz++)
//    {
//     //dprint(cz);
//     // index of interface to right of cell
//     const int ix = cx + 1;
//     const int iy = cy + 1;
//     const int iz = cz + 1;
//     {
//      // reference the 8 nodes to which we will
//      // write moment data for particles in this mesh cell.
//      //
//      arr1_double_fetch momentsArray[8];
//      arr2_double_fetch moments00 = moments[ix][iy];
//      arr2_double_fetch moments01 = moments[ix][cy];
//      arr2_double_fetch moments10 = moments[cx][iy];
//      arr2_double_fetch moments11 = moments[cx][cy];
//      momentsArray[0] = moments00[iz]; // moments000 
//      momentsArray[1] = moments00[cz]; // moments001 
//      momentsArray[2] = moments01[iz]; // moments010 
//      momentsArray[3] = moments01[cz]; // moments011 
//      momentsArray[4] = moments10[iz]; // moments100 
//      momentsArray[5] = moments10[cz]; // moments101 
//      momentsArray[6] = moments11[iz]; // moments110 
//      momentsArray[7] = moments11[cz]; // moments111 
//
//      // accumulator for moments per each of 8 threads
//      double momentsAcc[8][10][8];
//      const int numpcls_in_cell = pcls.get_numpcls_in_bucket(cx,cy,cz);
//      const int bucket_offset = pcls.get_bucket_offset(cx,cy,cz);
//      const int bucket_end = bucket_offset+numpcls_in_cell;
//
//      // data is not stride-1, so we do *not* use
//      // #pragma simd
//      {
//        // accumulators for moments per each of 8 threads
//        double momentsAcc[8][10];
//        memset(momentsAcc,0,sizeof(double)*8*10);
//        for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
//        {
//          const SpeciesParticle* pcl = &pcls.get_pcl(pidx);
//          // This depends on the fact that the memory
//          // occupied by a particle coincides with
//          // the alignment interval (64 bytes)
//          ALIGNED(pcl);
//          double velmoments[10];
//          double weights[8];
//          // compute the quadratic moments of velocity
//          //
//          const double ui=pcl->get_u();
//          const double vi=pcl->get_v();
//          const double wi=pcl->get_w();
//          const double uui=ui*ui;
//          const double uvi=ui*vi;
//          const double uwi=ui*wi;
//          const double vvi=vi*vi;
//          const double vwi=vi*wi;
//          const double wwi=wi*wi;
//          //double velmoments[10];
//          velmoments[0] = 1.;
//          velmoments[1] = ui;
//          velmoments[2] = vi;
//          velmoments[3] = wi;
//          velmoments[4] = uui;
//          velmoments[5] = uvi;
//          velmoments[6] = uwi;
//          velmoments[7] = vvi;
//          velmoments[8] = vwi;
//          velmoments[9] = wwi;
//        
//          // compute the weights to distribute the moments
//          //
//          //double weights[8];
//          const double abs_xpos = pcl->get_x();
//          const double abs_ypos = pcl->get_y();
//          const double abs_zpos = pcl->get_z();
//          const double rel_xpos = abs_xpos - xstart;
//          const double rel_ypos = abs_ypos - ystart;
//          const double rel_zpos = abs_zpos - zstart;
//          const double cxm1_pos = rel_xpos * inv_dx;
//          const double cym1_pos = rel_ypos * inv_dy;
//          const double czm1_pos = rel_zpos * inv_dz;
//          //if(true)
//          //{
//          //  const int cx_inf = int(floor(cxm1_pos));
//          //  const int cy_inf = int(floor(cym1_pos));
//          //  const int cz_inf = int(floor(czm1_pos));
//          //  assert_eq(cx-1,cx_inf);
//          //  assert_eq(cy-1,cy_inf);
//          //  assert_eq(cz-1,cz_inf);
//          //}
//          // fraction of the distance from the right of the cell
//          const double w1x = cx - cxm1_pos;
//          const double w1y = cy - cym1_pos;
//          const double w1z = cz - czm1_pos;
//          // fraction of distance from the left
//          const double w0x = 1-w1x;
//          const double w0y = 1-w1y;
//          const double w0z = 1-w1z;
//          // we are calculating a charge moment.
//          const double qi=pcl->get_q();
//          const double weight0 = qi*w0x;
//          const double weight1 = qi*w1x;
//          const double weight00 = weight0*w0y;
//          const double weight01 = weight0*w1y;
//          const double weight10 = weight1*w0y;
//          const double weight11 = weight1*w1y;
//          weights[0] = weight00*w0z; // weight000
//          weights[1] = weight00*w1z; // weight001
//          weights[2] = weight01*w0z; // weight010
//          weights[3] = weight01*w1z; // weight011
//          weights[4] = weight10*w0z; // weight100
//          weights[5] = weight10*w1z; // weight101
//          weights[6] = weight11*w0z; // weight110
//          weights[7] = weight11*w1z; // weight111
//        
//          // add moments for this particle
//          {
//            // which is the superior order for the following loop?
//            for(int c=0; c<8; c++)
//            for(int m=0; m<10; m++)
//            {
//              momentsAcc[c][m] += velmoments[m]*weights[c];
//            }
//          }
//        }
//        for(int c=0; c<8; c++)
//        for(int m=0; m<10; m++)
//        {
//          momentsArray[c][m] += momentsAcc[c][m];
//        }
//      }
//     }
//    }
//    #pragma omp master
//    { timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION); }
//
//    // reduction
//    #pragma omp master
//    { timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION); }
//    {
//      #pragma omp for collapse(2)
//      for(int i=0;i<nxn;i++){
//      for(int j=0;j<nyn;j++){
//      for(int k=0;k<nzn;k++)
//      {
//        rhons[is][i][j][k] = invVOL*moments[i][j][k][0];
//        Jxs  [is][i][j][k] = invVOL*moments[i][j][k][1];
//        Jys  [is][i][j][k] = invVOL*moments[i][j][k][2];
//        Jzs  [is][i][j][k] = invVOL*moments[i][j][k][3];
//        pXXsn[is][i][j][k] = invVOL*moments[i][j][k][4];
//        pXYsn[is][i][j][k] = invVOL*moments[i][j][k][5];
//        pXZsn[is][i][j][k] = invVOL*moments[i][j][k][6];
//        pYYsn[is][i][j][k] = invVOL*moments[i][j][k][7];
//        pYZsn[is][i][j][k] = invVOL*moments[i][j][k][8];
//        pZZsn[is][i][j][k] = invVOL*moments[i][j][k][9];
//      }}}
//    }
//    #pragma omp master
//    { timeTasks_end_task(TimeTasks::MOMENT_REDUCTION); }
//    // uncomment this and remove the loop below
//    // when we change to use asynchronous communication.
//    // communicateGhostP2G(is);
//  }
//  }
//  for (int i = 0; i < ns; i++)
//  {
//    communicateGhostP2G(i);
//  }
//}

// === end of methods_to_accumulate_moments ===

// This method assumes mirror boundary conditions;
// we therefore need to double the density on the boundary
// nodes to incorporate the mirror particles from the mirror
// cell just outside the domain.
//
/*! adjust densities on boundaries that are not periodic */
void EMfields3D::adjustNonPeriodicDensities(int is)
{
  const VirtualTopology3D *vct = &get_vct();
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][1][i][k] *= 2;
        Jxs  [is][1][i][k] *= 2;
        Jys  [is][1][i][k] *= 2;
        Jzs  [is][1][i][k] *= 2;
        pXXsn[is][1][i][k] *= 2;
        pXYsn[is][1][i][k] *= 2;
        pXZsn[is][1][i][k] *= 2;
        pYYsn[is][1][i][k] *= 2;
        pYZsn[is][1][i][k] *= 2;
        pZZsn[is][1][i][k] *= 2;
    }
  }
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][i][1][k] *= 2;
        Jxs  [is][i][1][k] *= 2;
        Jys  [is][i][1][k] *= 2;
        Jzs  [is][i][1][k] *= 2;
        pXXsn[is][i][1][k] *= 2;
        pXYsn[is][i][1][k] *= 2;
        pXZsn[is][i][1][k] *= 2;
        pYYsn[is][i][1][k] *= 2;
        pYZsn[is][i][1][k] *= 2;
        pZZsn[is][i][1][k] *= 2;
    }
  }
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
    {
        rhons[is][i][j][1] *= 2;
        Jxs  [is][i][j][1] *= 2;
        Jys  [is][i][j][1] *= 2;
        Jzs  [is][i][j][1] *= 2;
        pXXsn[is][i][j][1] *= 2;
        pXYsn[is][i][j][1] *= 2;
        pXZsn[is][i][j][1] *= 2;
        pYYsn[is][i][j][1] *= 2;
        pYZsn[is][i][j][1] *= 2;
        pZZsn[is][i][j][1] *= 2;
    }
  }
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][nxn - 2][i][k] *= 2;
        Jxs  [is][nxn - 2][i][k] *= 2;
        Jys  [is][nxn - 2][i][k] *= 2;
        Jzs  [is][nxn - 2][i][k] *= 2;
        pXXsn[is][nxn - 2][i][k] *= 2;
        pXYsn[is][nxn - 2][i][k] *= 2;
        pXZsn[is][nxn - 2][i][k] *= 2;
        pYYsn[is][nxn - 2][i][k] *= 2;
        pYZsn[is][nxn - 2][i][k] *= 2;
        pZZsn[is][nxn - 2][i][k] *= 2;
    }
  }
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][i][nyn - 2][k] *= 2;
        Jxs  [is][i][nyn - 2][k] *= 2;
        Jys  [is][i][nyn - 2][k] *= 2;
        Jzs  [is][i][nyn - 2][k] *= 2;
        pXXsn[is][i][nyn - 2][k] *= 2;
        pXYsn[is][i][nyn - 2][k] *= 2;
        pXZsn[is][i][nyn - 2][k] *= 2;
        pYYsn[is][i][nyn - 2][k] *= 2;
        pYZsn[is][i][nyn - 2][k] *= 2;
        pZZsn[is][i][nyn - 2][k] *= 2;
    }
  }
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
    {
        rhons[is][i][j][nzn - 2] *= 2;
        Jxs  [is][i][j][nzn - 2] *= 2;
        Jys  [is][i][j][nzn - 2] *= 2;
        Jzs  [is][i][j][nzn - 2] *= 2;
        pXXsn[is][i][j][nzn - 2] *= 2;
        pXYsn[is][i][j][nzn - 2] *= 2;
        pXZsn[is][i][j][nzn - 2] *= 2;
        pYYsn[is][i][j][nzn - 2] *= 2;
        pYZsn[is][i][j][nzn - 2] *= 2;
        pZZsn[is][i][j][nzn - 2] *= 2;
    }
  }
}

/*! communicate ghost for grid -> Particles interpolation */
void SpeciesMoms::communicateGhostP2G(int is)
{
  // interpolate adding common nodes among processors
  timeTasks_set_communicating();

  const VirtualTopology3D *vct = &get_vct();

  double ***moment0 = rhons.fetch_arr4()[is];
  double ***moment1 = Jxs  .fetch_arr4()[is];
  double ***moment2 = Jys  .fetch_arr4()[is];
  double ***moment3 = Jzs  .fetch_arr4()[is];
  double ***moment4 = pXXsn.fetch_arr4()[is];
  double ***moment5 = pXYsn.fetch_arr4()[is];
  double ***moment6 = pXZsn.fetch_arr4()[is];
  double ***moment7 = pYYsn.fetch_arr4()[is];
  double ***moment8 = pYZsn.fetch_arr4()[is];
  double ***moment9 = pZZsn.fetch_arr4()[is];
  // add the values for the shared nodes
  //
  communicateInterp(nxn, nyn, nzn, moment0, vct);
  communicateInterp(nxn, nyn, nzn, moment1, vct);
  communicateInterp(nxn, nyn, nzn, moment2, vct);
  communicateInterp(nxn, nyn, nzn, moment3, vct);
  communicateInterp(nxn, nyn, nzn, moment4, vct);
  communicateInterp(nxn, nyn, nzn, moment5, vct);
  communicateInterp(nxn, nyn, nzn, moment6, vct);
  communicateInterp(nxn, nyn, nzn, moment7, vct);
  communicateInterp(nxn, nyn, nzn, moment8, vct);
  communicateInterp(nxn, nyn, nzn, moment9, vct);
  // calculate the correct densities on the boundaries
  adjustNonPeriodicDensities(is);

  // populate the ghost nodes
  //
  communicateNode_P(nxn, nyn, nzn, moment0, vct);
  communicateNode_P(nxn, nyn, nzn, moment1, vct);
  communicateNode_P(nxn, nyn, nzn, moment2, vct);
  communicateNode_P(nxn, nyn, nzn, moment3, vct);
  communicateNode_P(nxn, nyn, nzn, moment4, vct);
  communicateNode_P(nxn, nyn, nzn, moment5, vct);
  communicateNode_P(nxn, nyn, nzn, moment6, vct);
  communicateNode_P(nxn, nyn, nzn, moment7, vct);
  communicateNode_P(nxn, nyn, nzn, moment8, vct);
  communicateNode_P(nxn, nyn, nzn, moment9, vct);
}

void MImoments::setZero()
{
  eprintf("check usage");

  for (register int i = 0; i < nxn; i++)
  for (register int j = 0; j < nyn; j++)
  for (register int k = 0; k < nzn; k++)
  {
    //Jx  [i][j][k] = 0.0;
    //Jy  [i][j][k] = 0.0;
    //Jz  [i][j][k] = 0.0;
    Jxh [i][j][k] = 0.0;
    Jyh [i][j][k] = 0.0;
    Jzh [i][j][k] = 0.0;
    rhon[i][j][k] = 0.0;
  }
  //for (register int i = 0; i < nxc; i++)
  //for (register int j = 0; j < nyc; j++)
  //for (register int k = 0; k < nzc; k++)
  //{
  //  rhoc[i][j][k] = 0.0;
  //  rhoh[i][j][k] = 0.0;
  //}
}

void SpeciesMoms::setZeroSpeciesMoms(int is)
{
  // set primary moments to zero
  //
  #pragma omp for collapse(2)
  for (register int i = 0; i < nxn; i++)
  for (register int j = 0; j < nyn; j++)
  for (register int k = 0; k < nzn; k++)
  {
    rhons[is][i][j][k] = 0.0;
    Jxs  [is][i][j][k] = 0.0;
    Jys  [is][i][j][k] = 0.0;
    Jzs  [is][i][j][k] = 0.0;
    pXXsn[is][i][j][k] = 0.0;
    pXYsn[is][i][j][k] = 0.0;
    pXZsn[is][i][j][k] = 0.0;
    pYYsn[is][i][j][k] = 0.0;
    pYZsn[is][i][j][k] = 0.0;
    pZZsn[is][i][j][k] = 0.0;
  }
}

/* sum the charge density of different species on nodes */
void EMfields3D::sumOverSpecies()
{
  for (int is = 0; is < ns; is++)
  for (register int i = 0; i < nxn; i++)
  for (register int j = 0; j < nyn; j++)
  for (register int k = 0; k < nzn; k++)
    rhon[i][j][k] += rhons[is][i][j][k];
}

/*! sum current density for different species */
void SpeciesMoms::sumOverSpeciesJ()
{
  for (int is = 0; is < ns; is++)
  for (register int i = 0; i < nxn; i++)
  for (register int j = 0; j < nyn; j++)
  for (register int k = 0; k < nzn; k++)
  {
    Jx[i][j][k] += Jxs[is][i][j][k];
    Jy[i][j][k] += Jys[is][i][j][k];
    Jz[i][j][k] += Jzs[is][i][j][k];
  }
}

arr3_double SpeciesMoms::ret_rhocs(int is) { return ca.get_N2C_no_ghosts(is, rhons); }
arr3_double SpeciesMoms::ret_Jxsc(int is) { return ca.get_N2C_no_ghosts(is, Jxs); }
arr3_double SpeciesMoms::ret_Jysc(int is) { return ca.get_N2C_no_ghosts(is, Jys); }
arr3_double SpeciesMoms::ret_Jzsc(int is) { return ca.get_N2C_no_ghosts(is, Jzs); }

// === section: MImoments methods ===

MImoments::compute_from_speciesMoms(const SpeciesMoms& speciesMoms,
  const_arr3_double Bx, const_arr3_double By, const_arr3_double Bz)
{
  //setZero();
  // sum all over the species
  //miMoments.sumOverSpecies();

  // copy the charge densities from the primitive moments
  const_arr3_double prim_rhons = speciesMoms.get_rhons();
  for(int is=0;is<ns;is++)
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    rhons[is][i][j][k] = prim_rhons[is][i][j][k];
  }
  if(Parameters::use_correct_smoothing())
  {
    for(int is=0;is<setting.col().getNs();is++)
    for(int i=0;i<Parameters::get_num_smoothings();i++)
      get_grid().smooth(rhons[is], 1);
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

  speciesMoms.calculateJhat(Jxh, Jyh, Jzh, Bx, By, Bz);
}

