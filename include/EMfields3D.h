/*!************************************************************************* EMfields3D.h - ElectroMagnetic fields definition ------------------- begin : May 2008 copyright : KU Leuven developers : Stefano Markidis, Giovanni Lapenta ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H

#include <iostream>
#include <sstream>

#include <math.h>
#include <mpi.h>

#include "hdf5.h"

#include "Alloc.h"
#include "Basic.h"
#include "Grid.h"
#include "TransArraySpace3D.h"
#include "CG.h"
#include "GMRES.h"
#include "Collective.h"
#include "ComNodes3D.h"
#include "ComInterpNodes3D.h"
//#include "TimeTasks.h"
#include "asserts.h"
#include "BCStructure.h"

using std::cout;
using std::cerr;
using std::endl;

/*! Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver @date May 2008 @par Copyright: (C) 2008 KUL @author Stefano Markidis, Giovanni Lapenta. @version 3.0 */

class Particles3Dcomm;
class Moments;
class EMfields3D                // :public Field
{
  public:
    /*! constructor */
    EMfields3D(Collective * col, Grid * grid);
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    void init(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! init beam */
    void initBEAM(VirtualTopology3D * vct, Grid * grid, Collective *col, double x_center, double y_center, double z_center, double radius);
    /*! initialize GEM challenge */
    void initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initOriginalGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialize GEM challenge with dipole-like tail without perturbation */
    void initGEMDipoleLikeTailNoPert(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialize GEM challenge with no Perturbation */
    void initGEMnoPert(VirtualTopology3D * vct, Grid * grid, Collective *col);
#ifdef BATSRUS
    /*! initialize from BATSRUS */
    void initBATSRUS(VirtualTopology3D * vct, Grid * grid, Collective * col);
#endif
    /*! Random initial field */
    void initRandomField(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! Init Force Free (JxB=0) */
    void initForceFree(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialized with rotated magnetic field */
    void initEM_rotate(VirtualTopology3D * vct, Grid * grid, Collective *col, double B, double theta);
    /*! add a perturbattion to charge density */
    void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid);
    /*! add a perturbattion to the EM field */
    void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid);
    /*! Initialise a combination of magnetic dipoles */
    void initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col);

    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Image of Poisson Solver (for SOLVER) */
    void PoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Image of Maxwell Solver (for Solver) */
    void MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Impose a constant charge inside a spherical zone of the domain */
    void ConstantChargePlanet(Grid * grid, VirtualTopology3D * vct, double R, double x_center, double y_center, double z_center);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBC(Grid * grid, VirtualTopology3D * vct);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBCv2(Grid * grid, VirtualTopology3D * vct);
    /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
    void calculateB(Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! fix B on the boundary for gem challange */
    void fixBgem(Grid * grid, VirtualTopology3D * vct);
    /*! fix B on the boundary for gem challange */
    void fixBforcefree(Grid * grid, VirtualTopology3D * vct);

    /*! Calculate the three components of Pi(implicit pressure) cross image vector */
    void PIdot(doubleArr3& PIdotX, doubleArr3& PIdotY, doubleArr3& PIdotZ,
      doubleCar3& vectX, doubleCar3& vectY, doubleCar3& vectZ, int ns, Grid * grid);
    /*! Calculate the three components of mu (implicit permeattivity) cross image vector */
    void MUdot(doubleArr3& MUdotX, doubleArr3& MUdotY, doubleArr3& MUdotZ,
      doubleCar3& vectX, doubleCar3& vectY, doubleCar3& vectZ, Grid * grid);
    /*! Calculate rho hat, Jx hat, Jy hat, Jz hat */
    void calculateHatFunctions(Grid * grid, VirtualTopology3D * vct);


    /*! communicate ghost for densities and interp rho from node to center */
    void interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid);
    /*! set to 0 all the densities fields */
    void setZeroDensities();
    /*! Sum rhon over species */
    void sumOverSpecies(VirtualTopology3D * vct);
    /*! Sum current over different species */
    void sumOverSpeciesJ();
    /*! Smoothing after the interpolation* */
    void smooth(double value, doubleArr3& vector, int type, Grid * grid, VirtualTopology3D * vct);
    /*! SPECIES: Smoothing after the interpolation for species fields* */
    void smooth(double value, doubleArr4& vector, int is, int type, Grid * grid, VirtualTopology3D * vct);
    /*! smooth the electric field */
    void smoothE(double value, VirtualTopology3D * vct, Collective *col);

    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct);
    void sumMoments(const Particles3Dcomm& pcls, Grid * grid, VirtualTopology3D * vct);
    /*! add accumulated moments to the moments for a given species */
    void addToSpeciesMoments(const Moments & in, int is);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction X to current density field at node X,Y,Z */
    void addJx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Y to current density field at node X,Y,Z */
    void addJy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Z to current density field at node X,Y,Z */
    void addJz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! add an amount of pressure density - direction XX to current density field at node X,Y,Z */
    void addPxx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XY to current density field at node X,Y,Z */
    void addPxy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
    void addPxz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YY to current density field at node X,Y,Z */
    void addPyy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
    void addPyz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
    void addPzz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! adjust densities on boundaries that are not periodic */
    void adjustNonPeriodicDensities(int is, VirtualTopology3D * vct);


    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(doubleArr3& imageX, doubleArr3& imageY, doubleArr3& imageZ,
      doubleCar3& vectorX, doubleCar3& vectorY, doubleCar3& vectorZ,
      int dir, Grid * grid);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(
      doubleArr3& imageX, doubleArr3& imageY, doubleArr3& imageZ,
      doubleCar3& vectorX,
      doubleCar3& vectorY,
      doubleCar3& vectorZ,
      int dir, Grid * grid);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(doubleArr3& vectorX, doubleArr3& vectorY, doubleArr3& vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(doubleArr3& vectorX, doubleArr3& vectorY, doubleArr3& vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorRightX(double **susxx, double **susyx, double **suszx);
    void sustensorLeftX (double **susxx, double **susyx, double **suszx);
    void sustensorRightY(double **susxy, double **susyy, double **suszy);
    void sustensorLeftY (double **susxy, double **susyy, double **suszy);
    void sustensorRightZ(double **susxz, double **susyz, double **suszz);
    void sustensorLeftZ (double **susxz, double **susyz, double **suszz);

    /*** accessor methods ***/

    /*! get Potential array */
    doubleArr3 getPHI() {return PHI;}

    // field components defined on nodes
    //
    double getEx(int X, int Y, int Z) const { return Ex.get(X,Y,Z);}
    double getEy(int X, int Y, int Z) const { return Ey.get(X,Y,Z);}
    double getEz(int X, int Y, int Z) const { return Ez.get(X,Y,Z);}
    double getBx(int X, int Y, int Z) const { return Bxn.get(X,Y,Z);}
    double getBy(int X, int Y, int Z) const { return Byn.get(X,Y,Z);}
    double getBz(int X, int Y, int Z) const { return Bzn.get(X,Y,Z);}
    //
    doubleArr3 getEx() { return Ex; }
    doubleArr3 getEy() { return Ey; }
    doubleArr3 getEz() { return Ez; }
    doubleArr3 getBx() { return Bxn; }
    doubleArr3 getBy() { return Byn; }
    doubleArr3 getBz() { return Bzn; }

    // field components without ghost cells
    //
    void getExc(doubleArr3& arr, Grid3DCU *grid);
    void getEyc(doubleArr3& arr, Grid3DCU *grid);
    void getEzc(doubleArr3& arr, Grid3DCU *grid);
    void getBxc(doubleArr3& arr);
    void getByc(doubleArr3& arr);
    void getBzc(doubleArr3& arr);

    doubleArr3 getRHOc() { return rhoc; }
    doubleArr3 getRHOn() { return rhon; }
    double getRHOc(int X, int Y, int Z) const { return rhoc.get(X,Y,Z);}
    double getRHOn(int X, int Y, int Z) const { return rhon.get(X,Y,Z);}

    // densities per species:
    //
    double getRHOcs(int X,int Y,int Z,int is)const{return rhocs.get(is,X,Y,Z);}
    double getRHOns(int X,int Y,int Z,int is)const{return rhons.get(is,X,Y,Z);}
    doubleArr4 getRHOns(){return rhons;}
    /* density on cells without ghost cells */
    void getRHOcs(doubleArr3& arr, Grid3DCU *grid, int is);

    double getBx_ext(int X, int Y, int Z) const{return Bx_ext.get(X,Y,Z);}
    double getBy_ext(int X, int Y, int Z) const{return By_ext.get(X,Y,Z);}
    double getBz_ext(int X, int Y, int Z) const{return Bz_ext.get(X,Y,Z);}
    
    doubleArr3 getBx_ext() { return Bx_ext; }
    doubleArr3 getBy_ext() { return By_ext; }
    doubleArr3 getBz_ext() { return Bz_ext; }

    doubleArr4 getpXXsn() { return pXXsn; }
    doubleArr4 getpXYsn() { return pXYsn; }
    doubleArr4 getpXZsn() { return pXZsn; }
    doubleArr4 getpYYsn() { return pYYsn; }
    doubleArr4 getpYZsn() { return pYZsn; }
    doubleArr4 getpZZsn() { return pZZsn; }

    double getJx(int X, int Y, int Z) const { return Jx.get(X,Y,Z);}
    double getJy(int X, int Y, int Z) const { return Jy.get(X,Y,Z);}
    double getJz(int X, int Y, int Z) const { return Jz.get(X,Y,Z);}
    doubleArr3 getJx() { return Jx; }
    doubleArr3 getJy() { return Jy; }
    doubleArr3 getJz() { return Jz; }
    doubleArr4 getJxs() { return Jxs; }
    doubleArr4 getJys() { return Jys; }
    doubleArr4 getJzs() { return Jzs; }

    double getJxs(int X,int Y,int Z,int is)const{return Jxs.get(is,X,Y,Z);}
    double getJys(int X,int Y,int Z,int is)const{return Jys.get(is,X,Y,Z);}
    double getJzs(int X,int Y,int Z,int is)const{return Jzs.get(is,X,Y,Z);}

    /*** accessor that require computing ***/

    // get current for species in all cells except ghost
    //
    void getJxsc(doubleArr3& arr, Grid3DCU *grid, int is);
    void getJysc(doubleArr3& arr, Grid3DCU *grid, int is);
    void getJzsc(doubleArr3& arr, Grid3DCU *grid, int is);

    /*! get the electric field energy */
    double getEenergy();
    /*! get the magnetic field energy */
    double getBenergy();

    /*! fetch array for summing moments of thread i */
    Moments& fetch_momentsArray(int i){
      assert_le(0,i);
      assert_le(i,sizeMomentsArray);
      return *momentsArray[i];
    }

    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col);

    /* ********************************* // VARIABLES ********************************* */
  private:
    /*! light speed */
    double c;
    /* 4*PI for normalization */
    double FourPI;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    /*! delt = c*th*dt */
    double delt;
    /*! number of particles species */
    int ns;
    /*! GEM challenge parameters */
    double B0x, B0y, B0z, delta;
    /** Earth Model parameters */
    double B1x, B1y, B1z;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! Boundary electron speed */
    double ue0, ve0, we0;


    // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
    /*! number of cells - X direction, including + 2 (guard cells) */
    int nxc;
    /*! number of nodes - X direction, including + 2 extra nodes for guard cells */
    int nxn;
    /*! number of cell - Y direction, including + 2 (guard cells) */
    int nyc;
    /*! number of nodes - Y direction, including + 2 extra nodes for guard cells */
    int nyn;
    /*! number of cell - Z direction, including + 2 (guard cells) */
    int nzc;
    /*! number of nodes - Z direction, including + 2 extra nodes for guard cells */
    int nzn;
    /*! local grid boundaries coordinate */
    double xStart, xEnd, yStart, yEnd, zStart, zEnd;
    /*! grid spacing */
    double dx, dy, dz, invVOL;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;
    /** source center - X direction   */
    double x_center;
    /** source center - Y direction   */
    double y_center;
    /** source center - Z direction   */
    double z_center;
    /** Characteristic length */
    double L_square;

    /*! PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes */
    doubleArray3 PHI;

    // Electric field components defined on nodes
    //
    doubleArray3 Ex;
    doubleArray3 Ey;
    doubleArray3 Ez;

    // implicit electric field components defined on nodes
    //
    doubleArray3 Exth;
    doubleArray3 Eyth;
    doubleArray3 Ezth;

    // magnetic field components defined on central points between nodes
    //
    doubleArray3 Bxc;
    doubleArray3 Byc;
    doubleArray3 Bzc;

    // magnetic field components defined on nodes
    //
    doubleArray3 Bxn;
    doubleArray3 Byn;
    doubleArray3 Bzn;

    // *************************************
    // TEMPORARY ARRAY
    // ************************************
    /*!some temporary arrays (for calculate hat functions) */
    doubleArray3 tempXC;
    doubleArray3 tempYC;
    doubleArray3 tempZC;
    doubleArray3 tempXN;
    doubleArray3 tempYN;
    doubleArray3 tempZN;
    /*! other temporary arrays (in MaxwellSource) */
    doubleArray3 tempC;
    doubleArray3 tempX;
    doubleArray3 tempY;
    doubleArray3 tempZ;
    doubleArray3 temp2X;
    doubleArray3 temp2Y;
    doubleArray3 temp2Z;
    /*! and some for MaxwellImage */
    doubleArray3 imageX;
    doubleArray3 imageY;
    doubleArray3 imageZ;
    doubleArray3 Dx;
    doubleArray3 Dy;
    doubleArray3 Dz;
    doubleArray3 vectX;
    doubleArray3 vectY;
    doubleArray3 vectZ;
    doubleArray3 divC;
    /* temporary arrays for summing moments */
    int sizeMomentsArray;
    Moments **momentsArray;


    // *******************************************************************************
    // *********** SOURCES **
    // *******************************************************************************

    /*! Charge density, defined on central points of the cell */
    doubleArray3 rhoc;
    /*! Charge density, defined on nodes */
    doubleArray3 rhon;
    /*! Implicit charge density, defined on central points of the cell */
    doubleArray3 rhoh;
    /*! SPECIES: charge density for each species, defined on nodes */
    doubleArray4 rhons;
    /*! SPECIES: charge density for each species, defined on central points of the cell */
    doubleArray4 rhocs;

    // current density defined on nodes
    //
    doubleArray3 Jx;
    doubleArray3 Jy;
    doubleArray3 Jz;

    // implicit current density defined on nodes
    //
    doubleArray3 Jxh;
    doubleArray3 Jyh;
    doubleArray3 Jzh;

    // species-specific current densities defined on nodes
    //
    doubleArray4 Jxs;
    doubleArray4 Jys;
    doubleArray4 Jzs;

    // magnetic field components defined on nodes
    //
    doubleArray3   Bx_ext;
    doubleArray3   By_ext;
    doubleArray3   Bz_ext;

    // external current, defined on nodes
    doubleArray3   Jx_ext;
    doubleArray3   Jy_ext;
    doubleArray3   Jz_ext;

    // pressure tensor components, defined on nodes
    doubleArray4 pXXsn;
    doubleArray4 pXYsn;
    doubleArray4 pXZsn;
    doubleArray4 pYYsn;
    doubleArray4 pYZsn;
    doubleArray4 pZZsn;

    /*! Field Boundary Condition
      0 = Dirichlet Boundary Condition: specifies the
          value on the boundary of the domain
      1 = Neumann Boundary Condition: specifies the value of
          derivative on the boundary of the domain
      2 = Periodic boundary condition */

    // boundary conditions for electrostatic potential
    //
    int bcPHIfaceXright;
    int bcPHIfaceXleft;
    int bcPHIfaceYright;
    int bcPHIfaceYleft;
    int bcPHIfaceZright;
    int bcPHIfaceZleft;

    /*! Boundary condition for electric field 0 = perfect conductor 1 = magnetic mirror */
    //
    // boundary conditions for EM field
    //
    int bcEMfaceXright;
    int bcEMfaceXleft;
    int bcEMfaceYright;
    int bcEMfaceYleft;
    int bcEMfaceZright;
    int bcEMfaceZleft;


    /*! GEM Challenge background ion */
    double *rhoINIT;
    /*! Drift of the species */
    bool *DriftSpecies;

    /*! boolean for divergence cleaning */
    bool PoissonCorrection;
    /*! RESTART BOOLEAN */
    int restart1;
    /*! String with the directory for the restart file */
    string RestartDirName;
    /*! Case */
    string Case;

    /*! CG tolerance criterium for stopping iterations */
    double CGtol;
    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;

    // OpenBC implementation

    injInfoFields *injFieldsLeft, *injFieldsRight, *injFieldsTop, *injFieldsBottom, *injFieldsFront, *injFieldsRear;

    injInfoFields* get_InfoFieldsLeft();
    injInfoFields* get_InfoFieldsTop();
    injInfoFields* get_InfoFieldsBottom();
    injInfoFields* get_InfoFieldsFront();
    injInfoFields* get_InfoFieldsRear();
    injInfoFields* get_InfoFieldsRight();

    void BoundaryConditionsB(doubleArr3& vectorX, doubleArr3& vectorY, doubleArr3& vectorZ,
      int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct);
    void BoundaryConditionsE(doubleArr3& vectorX, doubleArr3& vectorY, doubleArr3& vectorZ,
      int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct);
    void BoundaryConditionsEImage(doubleArr3& imageX, doubleArr3& imageY, doubleArr3& imageZ,
      const doubleCar3& vectorX, const doubleCar3& vectorY, const doubleCar3& vectorZ,
      int nx, int ny, int nz, VirtualTopology3D *vct,Grid *grid);
};

inline void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        rhons[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of charge density to current density - direction X to current density field on the node */
inline void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jxs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Y to current density field on the node */
inline void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jys[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Z to current density field on the node */
inline void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jzs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XX to current density field on the node */
inline void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXXsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XY to current density field on the node */
inline void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XZ to current density field on the node */
inline void EMfields3D::addPxz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YY to current density field on the node */
inline void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YZ to current density field on the node */
inline void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction ZZ to current density field on the node */
inline void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pZZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}


typedef EMfields3D Field;

#endif
