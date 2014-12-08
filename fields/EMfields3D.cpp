
#include <mpi.h>
#include "ipichdf5.h"
#include "EMfields3D.h"
#include "Basic.h"
#include "ComNodes3D.h"
#include "ComInterpNodes3D.h"
#include "CollectiveIO.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "CG.h"
#include "GMRES.h"
#include "BCStructure.h"
#include "Particles3Dcomm.h"
#include "TimeTasks.h"
#include "Moments.h"
#include "Parameters.h"
#include "ompdefs.h"
#include "debug.h"
#include "string.h" // for memset
#include "mic_basics.h"
#include "mic_particles.h"
#include "ipicmath.h" // for roundup_to_multiple
#include "Alloc.h"
#include "asserts.h"
#ifndef NO_HDF5
#include "Restart3D.h"
#endif

#include <iostream>
//#include <sstream>
using std::cout;
using std::endl;
using namespace iPic3D;

/*! destructor*/
EMfields3D::~EMfields3D() {
  delete [] qom;
  delete [] rhoINIT;
  delete injFieldsLeft;
  delete injFieldsRight;
  delete injFieldsTop;
  delete injFieldsBottom;
  delete injFieldsFront;
  delete injFieldsRear;
  delete Jx_ext,
  delete Jy_ext,
  delete Jz_ext,
  for(int i=0;i<sizeMomentsArray;i++) { delete moments10Array[i]; }
  delete [] moments10Array;
}

/*! constructor */
//
// We rely on the following rule from the C++ standard, section 12.6.2.5:
//
//   nonstatic data members shall be initialized in the order
//   they were declared in the class definition
//
// in particular, nxc, nyc, nzc and nxn, nyn, nzn are assumed
// initialized when subsequently used.
//
EMfields3D::EMfields3D(const Setting& setting_, Imoments *iMoments_)
: _setting(setting_),
  _col(setting_.col()),
  _vct(setting_.vct()),
  _grid(setting_.grid()),
  iMoments(iMoments_),
  nxc(grid->getNXC()),
  nxn(grid->getNXN()),
  nyc(grid->getNYC()),
  nyn(grid->getNYN()),
  nzc(grid->getNZC()),
  nzn(grid->getNZN()),
  ns(col->getNs()),
  c(col->getC()),
  dt(col->getDt()),
  th(col->getTh()),
  delt (c*th*dt), // declared after these
  //
  // array allocation: nodes
  //
  //fieldForPcls  (nxn, nyn, nzn, 2*DFIELD_3or4),
  Ex   (nxn, nyn, nzn),
  Ey   (nxn, nyn, nzn),
  Ez   (nxn, nyn, nzn),
  Exth (nxn, nyn, nzn), // eliminate?
  Eyth (nxn, nyn, nzn), // eliminate?
  Ezth (nxn, nyn, nzn), // eliminate?
  Bxn  (nxn, nyn, nzn),
  Byn  (nxn, nyn, nzn),
  Bzn  (nxn, nyn, nzn),
  // make these pointers that are only allocated if needed
  Bx_ext(nxn, nyn, nzn),
  By_ext(nxn, nyn, nzn),
  Bz_ext(nxn, nyn, nzn),
  // needed in MUdot
  Bx_tot(nxn, nyn, nzn),
  By_tot(nxn, nyn, nzn),
  Bz_tot(nxn, nyn, nzn),
  // smoothed components for particles or coarse moments
  Bx_smooth(nxn, nyn, nzn),
  By_smooth(nxn, nyn, nzn),
  Bz_smooth(nxn, nyn, nzn),
  Ex_smooth(nxn, nyn, nzn),
  Ey_smooth(nxn, nyn, nzn),
  Ez_smooth(nxn, nyn, nzn),

  // array allocation: central points 
  //
  PHI  (nxc, nyc, nzc),
  Bxc  (nxc, nyc, nzc),
  Byc  (nxc, nyc, nzc),
  Bzc  (nxc, nyc, nzc),
  rhoc (nxc, nyc, nzc),
  rhoh (nxc, nyc, nzc),

  // temporary arrays
  //
  tempC  (nxc, nyc, nzc),
  tempX  (nxn, nyn, nzn),
  tempY  (nxn, nyn, nzn),
  tempZ  (nxn, nyn, nzn),
  temp2X (nxn, nyn, nzn),
  temp2Y (nxn, nyn, nzn),
  temp2Z (nxn, nyn, nzn),
  imageX (nxn, nyn, nzn),
  imageY (nxn, nyn, nzn),
  imageZ (nxn, nyn, nzn),
  Dx (nxn, nyn, nzn),
  Dy (nxn, nyn, nzn),
  Dz (nxn, nyn, nzn),
  vectX (nxn, nyn, nzn),
  vectY (nxn, nyn, nzn),
  vectZ (nxn, nyn, nzn),
  divC  (nxc, nyc, nzc),
  // J_ext should not be allocated unless used.
  Jx_ext(0),
  Jy_ext(0),
  Jz_ext(0)
{
  // External imposed fields
  //
  Bx_ext.setall(col->getB1x());
  By_ext.setall(col->getB1y());
  Bz_ext.setall(col->getB1z());
  //
  PoissonCorrection = false;
  if (col->getPoissonCorrection()=="yes") PoissonCorrection = true;
  CGtol = col->getCGtol();
  GMREStol = col->getGMREStol();
  qom = new double[ns];
  for (int i = 0; i < ns; i++)
    qom[i] = col->getQOM(i);
  // boundary conditions: PHI and EM fields
  bcPHIfaceXright = col->getBcPHIfaceXright();
  bcPHIfaceXleft  = col->getBcPHIfaceXleft();
  bcPHIfaceYright = col->getBcPHIfaceYright();
  bcPHIfaceYleft  = col->getBcPHIfaceYleft();
  bcPHIfaceZright = col->getBcPHIfaceZright();
  bcPHIfaceZleft  = col->getBcPHIfaceZleft();

  bcEMfaceXright = col->getBcEMfaceXright();
  bcEMfaceXleft = col->getBcEMfaceXleft();
  bcEMfaceYright = col->getBcEMfaceYright();
  bcEMfaceYleft = col->getBcEMfaceYleft();
  bcEMfaceZright = col->getBcEMfaceZright();
  bcEMfaceZleft = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  delta = col->getDelta();
  // get the density background for the gem Challange
  rhoINIT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i = 0; i < ns; i++) {
    rhoINIT[i] = col->getRHOinit(i);
    if ((fabs(col->getW0(i)) != 0) || (fabs(col->getU0(i)) != 0)) // GEM and LHDI
      DriftSpecies[i] = true;
    else
      DriftSpecies[i] = false;
  }
  /*! parameters for GEM challenge */
  FourPI = 16 * atan(1.0);

  // OpenBC
  injFieldsLeft   = new injInfoFields(nxn, nyn, nzn);
  injFieldsRight  = new injInfoFields(nxn, nyn, nzn);
  injFieldsTop    = new injInfoFields(nxn, nyn, nzn);
  injFieldsBottom = new injInfoFields(nxn, nyn, nzn);
  injFieldsFront  = new injInfoFields(nxn, nyn, nzn);
  injFieldsRear   = new injInfoFields(nxn, nyn, nzn);
}

// === Section: methods_that_do_not_really_belong_in_this_class ===

/*! set to 0 all the densities fields */
//void EMfields3D::setZeroDensities() {
//  setZeroDerivedMoments();
//  setZeroPrimaryMoments();
//}

// === Section: methods_to_solve_fields ===

/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(arr3_double vectPhys, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
  for (register int j = 1; j < ny - 1; j++)
  for (register int k = 1; k < nz - 1; k++)
    vectPhys[i][j][k] = *vectSolver++;
}
/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(arr3_double vectPhys1, arr3_double vectPhys2, arr3_double vectPhys3, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vectPhys1[i][j][k] = *vectSolver++;
        vectPhys2[i][j][k] = *vectSolver++;
        vectPhys3[i][j][k] = *vectSolver++;
      }
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(double *vectSolver, const arr3_double vectPhys, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++)
        *vectSolver++ = vectPhys.get(i,j,k);
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(double *vectSolver, const arr3_double vectPhys1, const arr3_double vectPhys2, const arr3_double vectPhys3, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        *vectSolver++ = vectPhys1.get(i,j,k);
        *vectSolver++ = vectPhys2.get(i,j,k);
        *vectSolver++ = vectPhys3.get(i,j,k);
      }
}
/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::calculateE(const Imoments& iMoments)
{
  // define fields to use at open boundaries
  // based on magnetic field and Ohm's law
  // (this also updates boundaries for the magnetic field,
  // which is updated based on the updated electric field).
  updateInfoFields();

  const Collective *col = &get_col();
  const VirtualTopology3D * vct = &get_vct();
  const Grid *grid = &get_grid();

  if (get_vct().getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;

  array3_double divE     (nxc, nyc, nzc);
  array3_double gradPHIX (nxn, nyn, nzn);
  array3_double gradPHIY (nxn, nyn, nzn);
  array3_double gradPHIZ (nxn, nyn, nzn);

  double *xkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 E components
  double *bkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 components

  double *xkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  double *bkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  // set to zero all the stuff 
  eqValue(0.0, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, bkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, divE, nxc, nyc, nzc);
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, gradPHIX, nxn, nyn, nzn);
  eqValue(0.0, gradPHIY, nxn, nyn, nzn);
  eqValue(0.0, gradPHIZ, nxn, nyn, nzn);
  // Adjust E calculating laplacian(PHI) = div(E) -4*PI*rho DIVERGENCE CLEANING
  if (PoissonCorrection) {
    if (get_vct().getCartesian_rank() == 0)
      cout << "*** DIVERGENCE CLEANING ***" << endl;
    grid->divN2C(divE, Ex, Ey, Ez);
    scale(tempC, rhoc, -FourPI, nxc, nyc, nzc);
    sum(divE, tempC, nxc, nyc, nzc);
    // move to krylov space 
    phys2solver(bkrylovPoisson, divE, nxc, nyc, nzc);
    // use conjugate gradient first
    if (!CG(xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 3000, CGtol, &::PoissonImage, this)) {
      if (get_vct().getCartesian_rank() == 0)
        cout << "CG not Converged. Trying with GMRes. Consider to increase the number of the CG iterations" << endl;
      eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
      GMRES(&::PoissonImage, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 20, 200, GMREStol, this);
    }
    solver2phys(PHI, xkrylovPoisson, nxc, nyc, nzc);
    communicateCenterBC(nxc, nyc, nzc, PHI, 2, 2, 2, 2, 2, 2, vct);
    // calculate the gradient
    grid->gradC2N(gradPHIX, gradPHIY, gradPHIZ, PHI);
    // sub
    sub(Ex, gradPHIX, nxn, nyn, nzn);
    sub(Ey, gradPHIY, nxn, nyn, nzn);
    sub(Ez, gradPHIZ, nxn, nyn, nzn);

  }                             // end of divergence cleaning 
  if (get_vct().getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER ***" << endl;
  // prepare the source 
  MaxwellSource(bkrylov);
  phys2solver(xkrylov, Ex, Ey, Ez, nxn, nyn, nzn);
  // solver
  const void* registered_data[2] = {this, &iMoments.get_rhons()};
  GMRES(&::MaxwellImage, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2),
    bkrylov, 20, 200, GMREStol, registered_data);
  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, xkrylov, nxn, nyn, nzn);

  // update the electric field based on the average electric field
  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);

  // The original smoothing:
  // - smooths the evolved electric field E (incorrect)
  // - avoids smoothing the field used to update the fields (correct)
  // - uses the smoothed electric field used to push particles
  //
  // E rather than Eth is used to push particles in order to
  // damp artificial Chernikov radiation generated due to
  // the fact that the implicit moment method artificially
  // slows the speed of light.
  //
  if(Parameters::use_original_smoothing())
  {
    for(int i=0;i<Parameters::get_num_smoothings();i++)
    {
      grid->smooth(Ex, 1, get_col().get_bcEx());
      grid->smooth(Ey, 1, get_col().get_bcEy());
      grid->smooth(Ez, 1, get_col().get_bcEz());
    }
    // copy E to E_smooth
    for(int i=0;i<nxn;i++)
    for(int j=0;j<nyn;j++)
    for(int k=0;k<nzn;k++)
    {
      Ex_smooth[i][j][k] = Ex[i][j][k];
      Ey_smooth[i][j][k] = Ey[i][j][k];
      Ez_smooth[i][j][k] = Ez[i][j][k];
    }
  }

  // communicate so the interpolation can have good values
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx, vct);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy, vct);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz, vct);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx, vct);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy, vct);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz, vct);

  // OpenBC
  BoundaryConditionsE(Exth, Eyth, Ezth, nxn, nyn, nzn);
  BoundaryConditionsE(Ex, Ey, Ez, nxn, nyn, nzn);

  // correct smoothing:
  // - avoids smoothing the field used in subsequent field evolution
  // - uses the smoothed field used to push particles
  // - uses Eth rather than E for particles.
  //
  // note that theta=1 (i.e. Eth=E) is recommended to damp
  // artificial Chernikov radiation generated due to the
  // fact that the implicit moment method artificially slows
  // the speed of light.
  //
  if(Parameters::use_correct_smoothing())
  {
    // copy E to E_smooth
    for(int i=0;i<nxn;i++)
    for(int j=0;j<nyn;j++)
    for(int k=0;k<nzn;k++)
    {
      Ex_smooth[i][j][k] = Exth[i][j][k];
      Ey_smooth[i][j][k] = Eyth[i][j][k];
      Ez_smooth[i][j][k] = Ezth[i][j][k];
    }
    for(int i=0;i<Parameters::get_num_smoothings();i++)
    {
      grid->smooth(Ex_smooth, 1, get_col().bcEx);
      grid->smooth(Ey_smooth, 1, get_col().bcEy);
      grid->smooth(Ez_smooth, 1, get_col().bcEz);
    }
  }

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;
  delete[]xkrylovPoisson;
  delete[]bkrylovPoisson;
}

/*! Calculate source term for Maxwell solver */
void EMfields3D::MaxwellSource(double *bkrylov)
{
  const Collective *col = &get_col();
  const VirtualTopology3D * vct = &get_vct();
  const Grid *grid = &get_grid();

  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, tempXN, nxn, nyn, nzn);
  eqValue(0.0, tempYN, nxn, nyn, nzn);
  eqValue(0.0, tempZN, nxn, nyn, nzn);
  eqValue(0.0, temp2X, nxn, nyn, nzn);
  eqValue(0.0, temp2Y, nxn, nyn, nzn);
  eqValue(0.0, temp2Z, nxn, nyn, nzn);
  // communicate
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx, vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy, vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz, vct);

  if (get_col().getCase()=="ForceFree") fixBforcefree();
  if (get_col().getCase()=="GEM")       fixBgem();
  if (get_col().getCase()=="GEMnoPert") fixBgem();

  // OpenBC:
  BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc);

  // prepare curl of B for known term of Maxwell solver: for the source term
  grid->curlC2N(tempXN, tempYN, tempZN, Bxc, Byc, Bzc);
  scale(temp2X, iMoments.get_Jxh(), -FourPI / c, nxn, nyn, nzn);
  scale(temp2Y, iMoments.get_Jyh(), -FourPI / c, nxn, nyn, nzn);
  scale(temp2Z, iMoments.get_Jzh(), -FourPI / c, nxn, nyn, nzn);

  sum(temp2X, tempXN, nxn, nyn, nzn);
  sum(temp2Y, tempYN, nxn, nyn, nzn);
  sum(temp2Z, tempZN, nxn, nyn, nzn);
  scale(temp2X, delt, nxn, nyn, nzn);
  scale(temp2Y, delt, nxn, nyn, nzn);
  scale(temp2Z, delt, nxn, nyn, nzn);

  const int BCs[6] = {2,2,2,2,2,2};
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, BCs, vct);
  grid->gradC2N(tempX, tempY, tempZ, rhoh);

  scale(tempX, -delt * delt * FourPI, nxn, nyn, nzn);
  scale(tempY, -delt * delt * FourPI, nxn, nyn, nzn);
  scale(tempZ, -delt * delt * FourPI, nxn, nyn, nzn);
  // sum E, past values
  sum(tempX, Ex, nxn, nyn, nzn);
  sum(tempY, Ey, nxn, nyn, nzn);
  sum(tempZ, Ez, nxn, nyn, nzn);
  // sum curl(B) + jhat part
  sum(tempX, temp2X, nxn, nyn, nzn);
  sum(tempY, temp2Y, nxn, nyn, nzn);
  sum(tempZ, temp2Z, nxn, nyn, nzn);

  // Boundary condition in the known term
  // boundary condition: Xleft
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 0);
  // boundary condition: Xright
  if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 0);
  // boundary condition: Yleft
  if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 1);
  // boundary condition: Yright
  if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 1);
  // boundary condition: Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 2);
  // boundary condition: Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 2);

  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, nxn, nyn, nzn);
}

/*! Mapping of Maxwell image to give to solver */
//
// In the field solver, there is one layer of ghost cells.  The
// nodes on the ghost cells define two outer layers of nodes: the
// outermost nodes are clearly in the interior of the neighboring
// subdomain and can naturally be referred to as "ghost nodes",
// but the second-outermost layer is on the boundary between
// subdomains and thus does not clearly belong to any one process.
// Refer to these shared nodes as "boundary nodes".
//
// To compute the laplacian, we first compute the gradient
// at the center of each cell by differencing the values at
// the corners of the cell.  We then compute the laplacian
// (i.e. the divergence of the gradient) at each node by
// differencing the cell-center values in the cells sharing
// the node.
//
// The laplacian is required to be defined on all boundary
// and interior nodes.  
// 
// In the krylov solver, we make no attempt to use or to
// update the (outer) ghost nodes, and we assume (presumably
// correctly) that the boundary nodes are updated identically
// by all processes that share them.  Therefore, we must
// communicate gradient values in the ghost cells.  The
// subsequent computation of the divergence requires that
// this boundary communication first complete.
//
// An alternative way would be to communicate outer ghost node
// values after each update of Eth.  In this case, there would
// be no need for the 10=3*3+1 boundary communications in the body
// of MaxwellImage() entailed in the calls to lapN2N plus the
// call needed prior to the call to gradC2N.  Of course,
// we would then need to communicate the 3 components of the
// electric field for the outer ghost nodes prior to each call
// to MaxwellImage().  This second alternative would thus reduce
// communication by over a factor of 3.  Essentially, we would
// replace the cost of communicating cell-centered differences
// for ghost cell values with the cost of directly computing them.
//
// Also, while this second method does not increase the potential
// to avoid exposing latency, it can make it easier to do so.
//
// Another change that I would propose: define:
//
//   array4_double physical_vector(3,nxn,nyn,nzn);
//   arr3_double vectX = physical_vector[0];
//   arr3_double vectY = physical_vector[1];
//   arr3_double vectZ = physical_vector[2];
//   vector = &physical_vector[0][0][0][0];
//
// It is currently the case that boundary nodes are
// duplicated in "vector" and therefore receive a weight
// that is twice, four times, or eight times as much as
// other nodes in the Krylov inner product.  The definitions
// above would imply that ghost nodes also appear in the
// inner product.  To avoid this issue, we could simply zero
// ghost nodes before returning to the Krylov solver.  With
// the definitions above, phys2solver() would simply zero
// the ghost nodes and solver2phys() would populate them via
// communication.  Note that it would also be possible, if
// desired, to give duplicated nodes equal weight by
// rescaling their values in these two methods.
// 
// callback for generic GMRES solver
void ::MaxwellImage(double *im, double* vector, void** registered_data)
{
  EMfields3D* EMf = (EMfields3D*) registered_data[0];
  const vector_array3_double& rhons =
    *(const vector_array3_double*)(registered_data[1]);
  EMf->MaxwellImage(im, vector, rhons);
}
void EMfields3D::MaxwellImage(double *im, double* vector, const vector_array3_double& rhons)
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  eqValue(0.0, im, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, imageX, nxn, nyn, nzn);
  eqValue(0.0, imageY, nxn, nyn, nzn);
  eqValue(0.0, imageZ, nxn, nyn, nzn);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, Dx, nxn, nyn, nzn);
  eqValue(0.0, Dy, nxn, nyn, nzn);
  eqValue(0.0, Dz, nxn, nyn, nzn);
  // move from krylov space to physical space
  solver2phys(vectX, vectY, vectZ, vector, nxn, nyn, nzn);
  grid->lapN2N(imageX, vectX);
  grid->lapN2N(imageY, vectY);
  grid->lapN2N(imageZ, vectZ);
  neg(imageX, nxn, nyn, nzn);
  neg(imageY, nxn, nyn, nzn);
  neg(imageZ, nxn, nyn, nzn);
  // grad(div(mu dot E(n + theta)) mu dot E(n + theta) = D
  MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, rhons);
  grid->divN2C(divC, Dx, Dy, Dz);
  // communicate you should put BC 
  // think about the Physics 
  //const int BCs[6] = {1,1,1,1,1,1};
  const int BCs[6] = {2,2,2,2,2,2};
  // GO with Neumann, now then go with rho
  communicateCenterBC(nxc, nyc, nzc, divC, BCs, vct);

  grid->gradC2N(tempX, tempY, tempZ, divC);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta))
  sub(imageX, tempX, nxn, nyn, nzn);
  sub(imageY, tempY, nxn, nyn, nzn);
  sub(imageZ, tempZ, nxn, nyn, nzn);

  // scale delt*delt
  scale(imageX, delt * delt, nxn, nyn, nzn);
  scale(imageY, delt * delt, nxn, nyn, nzn);
  scale(imageZ, delt * delt, nxn, nyn, nzn);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta)) + eps dot E(n + theta)
  sum(imageX, Dx, nxn, nyn, nzn);
  sum(imageY, Dy, nxn, nyn, nzn);
  sum(imageZ, Dz, nxn, nyn, nzn);
  sum(imageX, vectX, nxn, nyn, nzn);
  sum(imageY, vectY, nxn, nyn, nzn);
  sum(imageZ, vectZ, nxn, nyn, nzn);

  // boundary condition: Xleft
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 0);
  // boundary condition: Xright
  if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 0);
  // boundary condition: Yleft
  if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 1);
  // boundary condition: Yright
  if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 1);
  // boundary condition: Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 2);
  // boundary condition: Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 2);

  // OpenBC
  BoundaryConditionsEImage(imageX, imageY, imageZ, vectX, vectY, vectZ, nxn, nyn, nzn);

  // move from physical space to krylov space
  phys2solver(im, imageX, imageY, imageZ, nxn, nyn, nzn);
}

/*! Calculate MU dot (vectX, vectY, vectZ) */
void EMfields3D::MUdot(
  arr3_double MUdotX,
  arr3_double MUdotY,
  arr3_double MUdotZ,
  const_arr3_double vectX,
  const_arr3_double vectY,
  const_arr3_double vectZ,
  const vector_array3_double& rhons)
{
  for (int i = 1; i < nxn - 1; i++)
  for (int j = 1; j < nyn - 1; j++)
  for (int k = 1; k < nzn - 1; k++)
  {
    MUdotX[i][j][k] = 0.0;
    MUdotY[i][j][k] = 0.0;
    MUdotZ[i][j][k] = 0.0;
  }
  for (int is = 0; is < ns; is++)
  {
    const double beta = .5 * qom[is] * dt / c;
    for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
    for (int k = 1; k < nzn - 1; k++)
    {
      const double omcx = beta * Bx_tot[i][j][k];
      const double omcy = beta * By_tot[i][j][k];
      const double omcz = beta * Bz_tot[i][j][k];
      const double edotb = vectX.get(i,j,k) * omcx + vectY.get(i,j,k) * omcy + vectZ.get(i,j,k) * omcz;
      const double denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
      MUdotX.fetch(i,j,k) += (vectX.get(i,j,k) + (vectY.get(i,j,k) * omcz - vectZ.get(i,j,k) * omcy + edotb * omcx)) * denom;
      MUdotY.fetch(i,j,k) += (vectY.get(i,j,k) + (vectZ.get(i,j,k) * omcx - vectX.get(i,j,k) * omcz + edotb * omcy)) * denom;
      MUdotZ.fetch(i,j,k) += (vectZ.get(i,j,k) + (vectX.get(i,j,k) * omcy - vectY.get(i,j,k) * omcx + edotb * omcz)) * denom;
    }
  }
}

/*! fix the B boundary when running gem */
void EMfields3D::fixBgem()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Bxc[i][nyc - 2][k] = Bxc[i][nyc - 1][k];
        Bxc[i][nyc - 3][k] = Bxc[i][nyc - 1][k];
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z;
        Bzc[i][nyc - 2][k] = B0z;
        Bzc[i][nyc - 3][k] = B0z;
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bxc[i][1][k] = Bxc[i][0][k];
        Bxc[i][2][k] = Bxc[i][0][k];
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z;
        Bzc[i][1][k] = B0z;
        Bzc[i][2][k] = B0z;
      }
  }
}

/*! fix the B boundary when running forcefree */
void EMfields3D::fixBforcefree()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z / cosh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 2][k] = B0z / cosh((grid->getYC(i, nyc - 2, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 3][k] = B0z / cosh((grid->getYC(i, nyc - 3, k) - Ly / 2) / delta);
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z / cosh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bzc[i][1][k] = B0z / cosh((grid->getYC(i, 1, k) - Ly / 2) / delta);
        Bzc[i][2][k] = B0z / cosh((grid->getYC(i, 2, k) - Ly / 2) / delta);
      }
  }
}

// === end of methods_to_solve_fields ===

// === Section: field_solver_routines ===

// update B_tot and B_smooth based on Bn and B_ext
//
void EMfield::update_total_B()
{
  // compute magnetic field used to push partices
  // or compute implicit moments
  {
    for(int i=0;i<nxn;i++)
    for(int j=0;j<nyn;j++)
    for(int k=0;k<nzn;k++)
    {
      // unsmoothed version: apply only to smoothed moments
      fetch_Bx_tot()[i][j][k] = Bxn[i][j][k]+fetch_Bx_ext()[i][j][k];
      fetch_By_tot()[i][j][k] = Byn[i][j][k]+fetch_By_ext()[i][j][k];
      fetch_Bz_tot()[i][j][k] = Bzn[i][j][k]+fetch_Bz_ext()[i][j][k];

      // smoothed version: apply only to unsmoothed moments or particles
      Bx_smooth[i][j][k] = fetch_Bx_tot()[i][j][k];
      By_smooth[i][j][k] = fetch_By_tot()[i][j][k];
      Bz_smooth[i][j][k] = fetch_Bz_tot()[i][j][k];
    }
    if(Parameters::use_correct_smoothing())
    {
      for(int i=0;i<Parameters::get_num_smoothings();i++)
      {
        grid->smooth(Bx_smooth, get_col().get_bcBx());
        grid->smooth(By_smooth, get_col().get_bcBy());
        grid->smooth(Bz_smooth, get_col().get_bcBz());
      }
    }
  }
}

/*! update the magnetic field using E(n+ theta) in Faraday's law */
void EMfields3D::advanceB()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getCartesian_rank() == 0)
    cout << "*** B CALCULATION ***" << endl;

  // calculate the curl of Eth
  grid->curlN2C(tempXC, tempYC, tempZC, Exth, Eyth, Ezth);
  // update the magnetic field
  addscale(-c * dt, 1, Bxc, tempXC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Byc, tempYC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Bzc, tempZC, nxc, nyc, nzc);
  // communicate ghost 
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx, vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy, vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz, vct);

  if (get_col().getCase()=="ForceFree") fixBforcefree();
  if (get_col().getCase()=="GEM")       fixBgem();
  if (get_col().getCase()=="GEMnoPert") fixBgem();

  // OpenBC:
  BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc);

  // interpolate C2N
  grid->interpC2N(Bxn, Bxc);
  grid->interpC2N(Byn, Byc);
  grid->interpC2N(Bzn, Bzc);

  communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx, vct);
  communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy, vct);
  communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz, vct);

  // update B_tot and B_smooth
  // (needed for particles, implicit moments, and MUdot).
  update_total_B();
}

// this has been split into calculateJhat and calculateRhoHat()
void EMfields3D::calculateHatFunctions()
{
  eprintf("defunct");
}
/*! Calculate rho hat */
void EMfields3D::calculateRhoHat(const Imoments& iMoments)
{
  const vector_array3_double& rhons = iMoments.get_rhons();
  const_arr3_double Jxh = iMoments.get_Jxh();
  const_arr3_double Jyh = iMoments.get_Jyh();
  const_arr3_double Jzh = iMoments.get_Jzh();
  // sum charge density over species and interpolate to cell centers
  // (also needed for Poisson solve)
  {
    array3_double rhon(nxn,nyn,nzn);
    rhon.setall(0.);
    for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
    for (register int k = 0; k < nzn; k++)
      rhon[i][j][k] += rhons[is][i][j][k];
    // calculate densities on centers from nodes
    get_grid().interpN2C(rhoc, rhon);
  }

  if(Parameters::use_original_smoothing())
  // otherwise rhon is already smoothed
  {
    smooth(rhoc, 0);
  }

  // calculate rho hat = rho - (dt*theta)div(jhat)
  get_grid().divN2C(tempXC, Jxh, Jyh, Jzh);
  scale(tempXC, -dt * th, nxc, nyc, nzc);
  sum(tempXC, rhoc, nxc, nyc, nzc);
  eq(rhoh, tempXC, nxc, nyc, nzc);
  // communicate rhoh
  const int BCs[6] = {2,2,2,2,2,2};
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, BCs, &get_vct());
}

// generic callback for GMRES or CG solver
void PoissonImage(double *image, double *vector, void** registered_data)
{
  EMfields3D* EMf = (EMfields3D*) registered_data[0];
  EMf->PoissonImage(image,vector);
}
/*! Image of Poisson Solver */
void EMfields3D::PoissonImage(double *image, double *vector)
{
  const Grid *grid = &get_grid();

  // allocate 2 three dimensional service vectors
  array3_double temp(nxc, nyc, nzc);
  array3_double im(nxc, nyc, nzc);
  eqValue(0.0, image, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, nxc, nyc, nzc);
  // calculate the laplacian
  grid->lapC2Cpoisson(im, temp);
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);
}

// === Section: initialization_routines ===

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
void MIsolver::initDipole(EMfields3D& EMf, Particles3Dcomm* part)
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

  // initialize rhons and E
  for (int i=0; i < nxn; i++)
  for (int j=0; j < nyn; j++)
  for (int k=0; k < nzn; k++)
  {
    for (int is=0; is < ns; is++){
      rhons[is][i][j][k] = rhoINIT[is]/FourPI;
    }
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

  for (int is=0 ; is<ns; is++)
  {
    grid->interpN2C(rhocs,is,rhons);
    part[i].maxwellian(rhocs);
  }

  if (col->getRestart_status()!=0)
    init_from_restart(rhocs);  // use the fields from restart file
  }
}

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

// === Section boundary_condition_routines ===

/*! Calculate the susceptibility on the boundary leftX */
void EMfields3D::sustensorLeftX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[1][j][k];
        omcy = beta * Byn[1][j][k];
        omcz = beta * Bzn[1][j][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][1][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary rightX */
void EMfields3D::sustensorRightX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[nxn - 2][j][k];
        omcy = beta * Byn[nxn - 2][j][k];
        omcz = beta * Bzn[nxn - 2][j][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][nxn - 2][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[i][1][k];
        omcy = beta * Byn[i][1][k];
        omcz = beta * Bzn[i][1][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][1][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[i][nyn - 2][k];
        omcy = beta * Byn[i][nyn - 2][k];
        omcz = beta * Bzn[i][nyn - 2][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][nyn - 2][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * Bxn[i][j][1];
        omcy = beta * Byn[i][j][1];
        omcz = beta * Bzn[i][j][1];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][1] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * Bxn[i][j][nzn - 2];
        omcy = beta * Byn[i][j][nzn - 2];
        omcz = beta * Bzn[i][j][nzn - 2];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][nyn - 2] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }
}

/*! Perfect conductor boundary conditions: LEFT wall */
void EMfields3D::perfectConductorLeft(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
  int dir)
{
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0:  // boundary condition on X-DIRECTION 
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorLeftX(susxx, susyx, suszx);
      for (int i=1; i <  nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[1][i][j] = vectorX.get(1,i,j) - (Ex[1][i][j] - susyx[i][j]*vectorY.get(1,i,j) - suszx[i][j]*vectorZ.get(1,i,j) - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[1][i][j] = vectorY.get(1,i,j) - 0.0*vectorY.get(2,i,j);
          imageZ[1][i][j] = vectorZ.get(1,i,j) - 0.0*vectorZ.get(2,i,j);
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorLeftY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[i][1][j] = vectorX.get(i,1,j) - 0.0*vectorX.get(i,2,j);
          imageY[i][1][j] = vectorY.get(i,1,j) - (Ey[i][1][j] - susxy[i][j]*vectorX.get(i,1,j) - suszy[i][j]*vectorZ.get(i,1,j) - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][1][j] = vectorZ.get(i,1,j) - 0.0*vectorZ.get(i,2,j);
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorLeftZ(susxy, susyy, suszy);
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          imageX[i][j][1] = vectorX.get(i,j,1);
          imageY[i][j][1] = vectorX.get(i,j,1);
          imageZ[i][j][1] = vectorZ.get(i,j,1) - (Ez[i][j][1] - susxz[i][j]*vectorX.get(i,j,1) - susyz[i][j]*vectorY.get(i,j,1) - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions: RIGHT wall */
void EMfields3D::perfectConductorRight(
  arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX,
  const_arr3_double vectorY,
  const_arr3_double vectorZ,
  int dir)
{
  double beta, omcx, omcy, omcz, denom;
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorRightX(susxx, susyx, suszx);
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[nxn-2][i][j] = vectorX.get(nxn-2,i,j) - (Ex[nxn-2][i][j] - susyx[i][j]*vectorY.get(nxn-2,i,j) - suszx[i][j]*vectorZ.get(nxn-2,i,j) - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[nxn-2][i][j] = vectorY.get(nxn-2,i,j) - 0.0 * vectorY.get(nxn-3,i,j);
          imageZ[nxn-2][i][j] = vectorZ.get(nxn-2,i,j) - 0.0 * vectorZ.get(nxn-3,i,j);
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);       
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorRightY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          imageX[i][nyn-2][j] = vectorX.get(i,nyn-2,j) - 0.0*vectorX.get(i,nyn-3,j);
          imageY[i][nyn-2][j] = vectorY.get(i,nyn-2,j) - (Ey[i][nyn-2][j] - susxy[i][j]*vectorX.get(i,nyn-2,j) - suszy[i][j]*vectorZ.get(i,nyn-2,j) - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][nyn-2][j] = vectorZ.get(i,nyn-2,j) - 0.0*vectorZ.get(i,nyn-3,j);
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION RIGHT
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorRightZ(susxz, susyz, suszz);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nyn-1;j++){
          imageX[i][j][nzn-2] = vectorX.get(i,j,nzn-2);
          imageY[i][j][nzn-2] = vectorY.get(i,j,nzn-2);
          imageZ[i][j][nzn-2] = vectorZ.get(i,j,nzn-2) - (Ez[i][j][nzn-2] - susxz[i][j]*vectorX.get(i,j,nzn-2) - susyz[i][j]*vectorY.get(i,j,nzn-2) - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);       
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::perfectConductorLeftS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  const double ue0 = col->getU0(0);
  const double ve0 = col->getV0(0);
  const double we0 = col->getW0(0);
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION LEFT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[1][i][j] = 0.0;
          vectorY[1][i][j] = ebc[1];
          vectorZ[1][i][j] = ebc[2];
          //+//          vectorX[1][i][j] = 0.0;
          //+//          vectorY[1][i][j] = 0.0;
          //+//          vectorZ[1][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][1][j] = ebc[0];
          vectorY[i][1][j] = 0.0;
          vectorZ[i][1][j] = ebc[2];
          //+//          vectorX[i][1][j] = 0.0;
          //+//          vectorY[i][1][j] = 0.0;
          //+//          vectorZ[i][1][j] = 0.0;
        }
      break;
    case 2: // boundary condition on Z-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][1] = ebc[0];
          vectorY[i][j][1] = ebc[1];
          vectorZ[i][j][1] = 0.0;
          //+//          vectorX[i][j][1] = 0.0;
          //+//          vectorY[i][j][1] = 0.0;
          //+//          vectorZ[i][j][1] = 0.0;
        }
      break;
  }
}

/*! Perfect conductor boundary conditions for source: RIGHT WALL */
void EMfields3D::perfectConductorRightS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  const double ue0 = col->getU0(0);
  const double ve0 = col->getV0(0);
  const double we0 = col->getW0(0);
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[nxn-2][i][j] = 0.0;
          vectorY[nxn-2][i][j] = ebc[1];
          vectorZ[nxn-2][i][j] = ebc[2];
          //+//          vectorX[nxn-2][i][j] = 0.0;
          //+//          vectorY[nxn-2][i][j] = 0.0;
          //+//          vectorZ[nxn-2][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][nyn-2][j] = ebc[0];
          vectorY[i][nyn-2][j] = 0.0;
          vectorZ[i][nyn-2][j] = ebc[2];
          //+//          vectorX[i][nyn-2][j] = 0.0;
          //+//          vectorY[i][nyn-2][j] = 0.0;
          //+//          vectorZ[i][nyn-2][j] = 0.0;
        }
      break;
    case 2:
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][nzn-2] = ebc[0];
          vectorY[i][j][nzn-2] = ebc[1];
          vectorZ[i][j][nzn-2] = 0.0;
          //+//          vectorX[i][j][nzn-2] = 0.0;
          //+//          vectorY[i][j][nzn-2] = 0.0;
          //+//          vectorZ[i][j][nzn-2] = 0.0;
        }
      break;
  }
}

// Open Boundary conditions implementation

void EMfields3D::updateInfoFields()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();

  double u_0, v_0, w_0;
  u_0=col->getU0(0);
  v_0=col->getV0(0);
  w_0=col->getW0(0);

  if (vct->getXleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< 3;i++)
      for (int j=0; j<nyn;j++)
        for (int k=0; k<nzn;k++){

          injFieldsLeft->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsLeft->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsLeft->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsLeft->BxITemp[i][j][k]=B0x;
          injFieldsLeft->ByITemp[i][j][k]=B0y;
          injFieldsLeft->BzITemp[i][j][k]=B0z;
        }
  }

  if (vct->getXright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=nxn-3; i< nxn; i++)
      for (int j=0; j<nyn; j++)
        for (int k=0; k<nzn; k++){

          injFieldsRight->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsRight->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsRight->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsRight->BxITemp[i][j][k]=B0x;
          injFieldsRight->ByITemp[i][j][k]=B0y;
          injFieldsRight->BzITemp[i][j][k]=B0z;

        }

  }

  if (vct->getYleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<3;j++)
        for (int k=0; k<nzn;k++){

          injFieldsBottom->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsBottom->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsBottom->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsBottom->BxITemp[i][j][k]=B0x;
          injFieldsBottom->ByITemp[i][j][k]=B0y;
          injFieldsBottom->BzITemp[i][j][k]=B0z;
        }

  }
  if (vct->getYright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=nyn-3; j<nyn;j++)
        for (int k=0; k<nzn;k++){

          injFieldsTop->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsTop->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsTop->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsTop->BxITemp[i][j][k]=B0x;
          injFieldsTop->ByITemp[i][j][k]=B0y;
          injFieldsTop->BzITemp[i][j][k]=B0z;
        }

  }
  if (vct->getZleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<nyn;j++)
        for (int k=0; k<3;k++){

          injFieldsRear->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsRear->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsRear->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsRear->BxITemp[i][j][k]=B0x;
          injFieldsRear->ByITemp[i][j][k]=B0y;
          injFieldsRear->BzITemp[i][j][k]=B0z;
        }

  }

  if (vct->getZright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<nyn;j++)
        for (int k=nzn-3; k<nzn;k++){

          injFieldsFront->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsFront->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsFront->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsFront->BxITemp[i][j][k]=B0x;
          injFieldsFront->ByITemp[i][j][k]=B0y;
          injFieldsFront->BzITemp[i][j][k]=B0z;
        }
  }

}

void EMfields3D::BoundaryConditionsEImage(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[0][j][k] = vectorX[0][j][k] - injFieldsLeft->ExITemp[0][j][k];
        imageY[0][j][k] = vectorY[0][j][k] - injFieldsLeft->EyITemp[0][j][k];
        imageZ[0][j][k] = vectorZ[0][j][k] - injFieldsLeft->EzITemp[0][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[nx-1][j][k] = vectorX[nx-1][j][k]- injFieldsRight->ExITemp[nx-1][j][k];
        imageY[nx-1][j][k] = vectorY[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];
        imageZ[nx-1][j][k] = vectorZ[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];

      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][0][k] = vectorX[i][0][k]-injFieldsBottom->ExITemp[i][0][k];
        imageY[i][0][k] = vectorY[i][0][k]-injFieldsBottom->EyITemp[i][0][k];
        imageZ[i][0][k] = vectorZ[i][0][k]-injFieldsBottom->EzITemp[i][0][k];
      }

  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][ny-1][k] = vectorX[i][ny-1][k]-injFieldsTop->ExITemp[i][nx-1][k];
        imageY[i][ny-1][k] = vectorY[i][ny-1][k]-injFieldsTop->EyITemp[i][nx-1][k];
        imageZ[i][ny-1][k] = vectorZ[i][ny-1][k]-injFieldsTop->EzITemp[i][nx-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][0] = vectorX[i][j][0]-injFieldsFront->ExITemp[i][j][0];
        imageY[i][j][0] = vectorY[i][j][0]-injFieldsFront->EyITemp[i][j][0];
        imageZ[i][j][0] = vectorZ[i][j][0]-injFieldsFront->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][nz-1] = vectorX[i][j][nz-1]-injFieldsRear->ExITemp[i][j][nz-1];
        imageY[i][j][nz-1] = vectorY[i][j][nz-1]-injFieldsRear->EyITemp[i][j][nz-1];
        imageZ[i][j][nz-1] = vectorZ[i][j][nz-1]-injFieldsRear->EzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsB(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
        vectorX[0][j][k] = injFieldsLeft->BxITemp[0][j][k];
        vectorY[0][j][k] = injFieldsLeft->ByITemp[0][j][k];
        vectorZ[0][j][k] = injFieldsLeft->BzITemp[0][j][k];

//      vectorX[1][j][k] = injFieldsLeft->BxITemp[1][j][k];
//      vectorY[1][j][k] = injFieldsLeft->ByITemp[1][j][k];
//      vectorZ[1][j][k] = injFieldsLeft->BzITemp[1][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
//      vectorX[nx-2][j][k] = injFieldsRight->BxITemp[nx-2][j][k];
//      vectorY[nx-2][j][k] = injFieldsRight->ByITemp[nx-2][j][k];
//      vectorZ[nx-2][j][k] = injFieldsRight->BzITemp[nx-2][j][k];

        vectorX[nx-1][j][k] = injFieldsRight->BxITemp[nx-1][j][k];
        vectorY[nx-1][j][k] = injFieldsRight->ByITemp[nx-1][j][k];
        vectorZ[nx-1][j][k] = injFieldsRight->BzITemp[nx-1][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][1][k] = injFieldsBottom->BxITemp[i][1][k];
//      vectorY[i][1][k] = injFieldsBottom->ByITemp[i][1][k];
//      vectorZ[i][1][k] = injFieldsBottom->BzITemp[i][1][k];

        vectorX[i][0][k] = injFieldsBottom->BxITemp[i][0][k];
        vectorY[i][0][k] = injFieldsBottom->ByITemp[i][0][k];
        vectorZ[i][0][k] = injFieldsBottom->BzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][ny-2][k] = injFieldsTop->BxITemp[i][ny-2][k];
//      vectorY[i][ny-2][k] = injFieldsTop->ByITemp[i][ny-2][k];
//      vectorZ[i][ny-2][k] = injFieldsTop->BzITemp[i][ny-2][k];

        vectorX[i][ny-1][k] = injFieldsTop->BxITemp[i][ny-1][k];
        vectorY[i][ny-1][k] = injFieldsTop->ByITemp[i][ny-1][k];
        vectorZ[i][ny-1][k] = injFieldsTop->BzITemp[i][ny-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][1] = injFieldsRear->BxITemp[i][j][1];
//      vectorY[i][j][1] = injFieldsRear->ByITemp[i][j][1];
//      vectorZ[i][j][1] = injFieldsRear->BzITemp[i][j][1];

        vectorX[i][j][0] = injFieldsRear->BxITemp[i][j][0];
        vectorY[i][j][0] = injFieldsRear->ByITemp[i][j][0];
        vectorZ[i][j][0] = injFieldsRear->BzITemp[i][j][0];
      }
  }


  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][nz-2] = injFieldsFront->BxITemp[i][j][nz-2];
//      vectorY[i][j][nz-2] = injFieldsFront->ByITemp[i][j][nz-2];
//      vectorZ[i][j][nz-2] = injFieldsFront->BzITemp[i][j][nz-2];

        vectorX[i][j][nz-1] = injFieldsFront->BxITemp[i][j][nz-1];
        vectorY[i][j][nz-1] = injFieldsFront->ByITemp[i][j][nz-1];
        vectorZ[i][j][nz-1] = injFieldsFront->BzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsE(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
        vectorX[1][j][k] = injFieldsLeft->ExITemp[1][j][k];
        vectorY[1][j][k] = injFieldsLeft->EyITemp[1][j][k];
        vectorZ[1][j][k] = injFieldsLeft->EzITemp[1][j][k];

//      vectorX[0][j][k] = injFieldsLeft->ExITemp[0][j][k];
//      vectorY[0][j][k] = injFieldsLeft->EyITemp[0][j][k];
//      vectorZ[0][j][k] = injFieldsLeft->EzITemp[0][j][k];
      } 
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){

//      vectorX[nx-2][j][k] = injFieldsRight->ExITemp[1][j][k];
//      vectorY[nx-2][j][k] = injFieldsRight->EyITemp[1][j][k];
//      vectorZ[nx-2][j][k] = injFieldsRight->EzITemp[1][j][k];

        vectorX[nx-1][j][k] = injFieldsRight->ExITemp[0][j][k];
        vectorY[nx-1][j][k] = injFieldsRight->EyITemp[0][j][k];
        vectorZ[nx-1][j][k] = injFieldsRight->EzITemp[0][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][1][k] = injFieldsBottom->ExITemp[i][1][k];
//      vectorY[i][1][k] = injFieldsBottom->EyITemp[i][1][k];
//      vectorZ[i][1][k] = injFieldsBottom->EzITemp[i][1][k];

        vectorX[i][0][k] = injFieldsBottom->ExITemp[i][0][k];
        vectorY[i][0][k] = injFieldsBottom->EyITemp[i][0][k];
        vectorZ[i][0][k] = injFieldsBottom->EzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][ny-2][k] = injFieldsTop->ExITemp[i][1][k];
//      vectorY[i][ny-2][k] = injFieldsTop->EyITemp[i][1][k];
//      vectorZ[i][ny-2][k] = injFieldsTop->EzITemp[i][1][k];

        vectorX[i][ny-1][k] = injFieldsTop->ExITemp[i][0][k];
        vectorY[i][ny-1][k] = injFieldsTop->EyITemp[i][0][k];
        vectorZ[i][ny-1][k] = injFieldsTop->EzITemp[i][0][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][1] = injFieldsRear->ExITemp[i][j][1];
//      vectorY[i][j][1] = injFieldsRear->EyITemp[i][j][1];
//      vectorZ[i][j][1] = injFieldsRear->EzITemp[i][j][1];

        vectorX[i][j][0] = injFieldsRear->ExITemp[i][j][0];
        vectorY[i][j][0] = injFieldsRear->EyITemp[i][j][0];
        vectorZ[i][j][0] = injFieldsRear->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][nz-2] = injFieldsFront->ExITemp[i][j][1];
//      vectorY[i][j][nz-2] = injFieldsFront->EyITemp[i][j][1];
//      vectorZ[i][j][nz-2] = injFieldsFront->EzITemp[i][j][1];

        vectorX[i][j][nz-1] = injFieldsFront->ExITemp[i][j][0];
        vectorY[i][j][nz-1] = injFieldsFront->EyITemp[i][j][0];
        vectorZ[i][j][nz-1] = injFieldsFront->EzITemp[i][j][0];
      }
  }
}

// === end boundary_condition_routines ===

// === Section: background_conditions_routines ===

void EMfields3D::ConstantChargeOpenBCv2()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

    ff = 1.0;
    if (is == 0) ff = -1.0;

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = rhons[is][4][j][k];
          rhons[is][1][j][k] = rhons[is][4][j][k];
          rhons[is][2][j][k] = rhons[is][4][j][k];
          rhons[is][3][j][k] = rhons[is][4][j][k];
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-3][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-2][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-1][j][k] = rhons[is][nx-5][j][k];
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = rhons[is][i][4][k];
          rhons[is][i][1][k] = rhons[is][i][4][k];
          rhons[is][i][2][k] = rhons[is][i][4][k];
          rhons[is][i][3][k] = rhons[is][i][4][k];
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-3][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-2][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-1][k] = rhons[is][i][ny-5][k];
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = rhons[is][i][j][4];
          rhons[is][i][j][1] = rhons[is][i][j][4];
          rhons[is][i][j][2] = rhons[is][i][j][4];
          rhons[is][i][j][3] = rhons[is][i][j][4];
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-3] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-2] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-1] = rhons[is][i][j][nz-5];
        }
    }
  }

}

void Imoments::ConstantChargeOpenBC()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

    ff = 1.0;
    if (is == 0) ff = -1.0;

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcEMfaceXleft ==2)) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][1][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][3][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcEMfaceXright ==2)) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-3][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-1][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcEMfaceYleft ==2))  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][1][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][3][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcEMfaceYright ==2))  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-3][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-1][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && (bcEMfaceZleft ==2))  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][1] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][3] = ff * rhoINIT[is] / FourPI;
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && (bcEMfaceZright ==2))  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-3] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-1] = ff * rhoINIT[is] / FourPI;
        }
    }
  }

}

void Imoments::ConstantChargePlanet()
{
  const Collective *col = &setting.get_col();
  const double R = col->getL_square();
  const double x_center = col->getx_center()
  const double y_center = col->gety_center()
  const double z_center = col->getz_center()

  const Grid *grid = &setting.get_grid();

  for (int is = 0; is < ns; is++)
  {
    double ff = 1.0;
    if (is == 0) ff = -1.0;
    for (int i = 1; i < nxn; i++)
    for (int j = 1; j < nyn; j++)
    for (int k = 1; k < nzn; k++)
    {
       const double xd = grid->getXN(i) - x_center;
       const double yd = grid->getYN(j) - y_center;
       const double zd = grid->getZN(k) - z_center;

       if ((xd*xd+yd*yd+zd*zd) <= R*R) {
         rhons[is][i][j][k] = ff * rhoINIT[is] / FourPI;
       }
    }
  }
}

// === end background_conditions_routines ===

// === Section output_reporting_routines ===

arr3_double EMfields3D::ret_Bxc() { return ca.get_no_ghosts(Bxc); }
arr3_double EMfields3D::ret_Byc() { return ca.get_no_ghosts(Byc); }
arr3_double EMfields3D::ret_Bzc() { return ca.get_no_ghosts(Bzc); }
arr3_double EMfields3D::ret_Exc() { return ca.get_N2C_no_ghosts(Ex); }
arr3_double EMfields3D::ret_Eyc() { return ca.get_N2C_no_ghosts(Ey); }
arr3_double EMfields3D::ret_Ezc() { return ca.get_N2C_no_ghosts(Ez); }

/*! get the electric field energy */
double EMfields3D::getEenergy(void) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++)
        localEenergy += .5 * grid->getVOL() * (Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI);

  MPI_Allreduce(&localEenergy, &totalEenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalEenergy);

}
/*! get the magnetic field energy */
double EMfields3D::getBenergy(void) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++){
        Bxt = Bxn[i][j][k]+get_Bx_ext(i,j,k);
        Byt = Byn[i][j][k]+get_By_ext(i,j,k);
        Bzt = Bzn[i][j][k]+get_Bz_ext(i,j,k);
        localBenergy += .5*grid->getVOL()*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalBenergy);
}

// === Section output_reporting_routines ===

