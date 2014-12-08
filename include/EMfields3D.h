/*!************************************************************************* EMfields3D.h - ElectroMagnetic fields definition ------------------- begin : May 2008 copyright : KU Leuven developers : Stefano Markidis, Giovanni Lapenta ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H

#include "asserts.h"
//#include "BCStructure.h"
#include "ipicfwd.h"
#include "Alloc.h"
struct injInfoFields;

/*! Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver @date May 2008 @par Copyright: (C) 2008 KUL @author Stefano Markidis, Giovanni Lapenta. @version 3.0 */

class EMfields3D
{
  public:
    /*! constructor */
    EMfields3D(Collective * col, Grid * grid, VirtualTopology3D *vct);
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    virtual void init();

    /*! init beam */
    void initBEAM(double x_center, double y_center, double z_center, double radius);
    /*! initialize GEM challenge */
    void initGEM();
    void initOriginalGEM();
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation();
    /*! initialize GEM challenge with dipole-like tail without perturbation */
    void initGEMDipoleLikeTailNoPert();
    /*! initialize GEM challenge with no Perturbation */
    void initGEMnoPert();
#ifdef BATSRUS
    /*! initialize from BATSRUS */
    void initBATSRUS();
#endif
    /*! Random initial field */
    void initRandomField();
    /*! Init Force Free (JxB=0) */
    void initForceFree();
    /*! initialized with rotated magnetic field */
    void initEM_rotate(double B, double theta);
    /*! add a perturbattion to charge density */
    void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid);
    /*! add a perturbattion to the EM field */
    void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid);
    /*! Initialise a combination of magnetic dipoles */
    void initDipole();

    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE();
    /*! Image of Poisson Solver (for SOLVER) */
    void PoissonImage(double *image, double *vector);
    /*! Image of Maxwell Solver (for Solver) */
    void MaxwellImage(double *im, double *vector);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov);
    /*! Impose a constant charge inside a spherical zone of the domain */
    void ConstantChargePlanet(double R, double x_center, double y_center, double z_center);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBC();
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBCv2();
    /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
    void advanceB();
    /*! fix B on the boundary for gem challange */
    void fixBgem();
    /*! fix B on the boundary for gem challange */
    void fixBforcefree();

    /*! Calculate the three components of mu (implicit permeattivity)
        cross image vector */
    // access to this is needed in the fields solver
    void MUdot(arr3_double MUdotX, arr3_double MUdotY, arr3_double MUdotZ,
      const_arr3_double vectX, const_arr3_double vectY, const_arr3_double vectZ);

    /*! copy the field data to the array used to move the particles */
    void set_fieldForPcls();

    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
      int dir);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(
      arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX,
      const_arr3_double vectorY,
      const_arr3_double vectorZ,
      int dir);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorRightX(double **susxx, double **susyx, double **suszx);
    void sustensorLeftX (double **susxx, double **susyx, double **suszx);
    void sustensorRightY(double **susxy, double **susyy, double **suszy);
    void sustensorLeftY (double **susxy, double **susyy, double **suszy);
    void sustensorRightZ(double **susxz, double **susyz, double **suszz);
    void sustensorLeftZ (double **susxz, double **susyz, double **suszz);

    /*** accessor methods ***/

    /*! get Potential array */
    arr3_double getPHI() {return PHI;}

    // field components defined on nodes
    //
    double getEx(int X, int Y, int Z) const { return Ex.get(X,Y,Z);}
    double getEy(int X, int Y, int Z) const { return Ey.get(X,Y,Z);}
    double getEz(int X, int Y, int Z) const { return Ez.get(X,Y,Z);}
    double getBx(int X, int Y, int Z) const { return Bxn.get(X,Y,Z);}
    double getBy(int X, int Y, int Z) const { return Byn.get(X,Y,Z);}
    double getBz(int X, int Y, int Z) const { return Bzn.get(X,Y,Z);}
    //
    //arr3_double getEx() { return Ex; }
    //arr3_double getEy() { return Ey; }
    //arr3_double getEz() { return Ez; }
    //arr3_double getBx() { return Bxn; }
    //arr3_double getBy() { return Byn; }
    //arr3_double getBz() { return Bzn; }

    arr3_double fetch_Ex() { return Ex; }
    arr3_double fetch_Ey() { return Ey; }
    arr3_double fetch_Ez() { return Ez; }
    arr3_double fetch_Bxn() { return Bxn; }
    arr3_double fetch_Byn() { return Byn; }
    arr3_double fetch_Bzn() { return Bzn; }
    arr3_double fetch_Bxc() { return Bxc; }
    arr3_double fetch_Byc() { return Byc; }
    arr3_double fetch_Bzc() { return Bzc; }
    arr3_double fetch_Bx_ext() { return Bx_ext; }
    arr3_double fetch_By_ext() { return By_ext; }
    arr3_double fetch_Bz_ext() { return Bz_ext; }

    const_arr3_double get_Ex()const{ return Ex; }
    const_arr3_double get_Ey()const{ return Ey; }
    const_arr3_double get_Ez()const{ return Ez; }
    const_arr3_double get_Bxn()const{ return Bxn; }
    const_arr3_double get_Byn()const{ return Byn; }
    const_arr3_double get_Bzn()const{ return Bzn; }
    const_arr3_double get_Bxc()const{ return Bxc; }
    const_arr3_double get_Byc()const{ return Byc; }
    const_arr3_double get_Bzc()const{ return Bzc; }
    const_arr3_double get_Bx_smooth()const{ return Bx_smooth; }
    const_arr3_double get_By_smooth()const{ return By_smooth; }
    const_arr3_double get_Bz_smooth()const{ return Bz_smooth; }

    // field components without ghost cells
    //
    //arr3_double getExc();
    //arr3_double getEyc();
    //arr3_double getEzc();
    //arr3_double getBxc();
    //arr3_double getByc();
    //arr3_double getBzc();

    /*! get the electric field energy */
    double getEenergy();
    /*! get the magnetic field energy */
    double getBenergy();

    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields();

  public: // accessors
    const Setting& get_setting()const{return _setting;}
    const Collective& get_col()const{return _col;}
    const Grid& get_grid()const{return _grid;};
    const VirtualTopology3D& get_vct()const{return _vct;}
    /* ********************************* // VARIABLES ********************************* */

  private:
    // access to global data
    const Setting& _setting;
    const Collective& _col;
    const VirtualTopology3D&_vct;
    const Grid& _grid;
    //
    // copies of Collective parameters
    // (We should try to get ride of these and instead
    // access this stuff via get_col().)
    //
    /*! light speed */
    double c;
    /* 4*PI for normalization */
    double FourPI;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
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

    /*! PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes */
    array3_double PHI;

    // Electric field components defined on nodes
    //
    array3_double Ex;
    array3_double Ey;
    array3_double Ez;

    // implicit electric field components defined on nodes
    //
    array3_double Exth;
    array3_double Eyth;
    array3_double Ezth;

    // magnetic field components defined on central points between nodes
    //
    array3_double Bxc;
    array3_double Byc;
    array3_double Bzc;

    // magnetic field components defined on nodes
    //
    array3_double Bxn;
    array3_double Byn;
    array3_double Bzn;
    array3_double* Bxtot;

    // *************************************
    // TEMPORARY ARRAY
    // ************************************
    /*! other temporary arrays (in MaxwellSource) */
    //array3_double tempC;
    //array3_double tempX;
    //array3_double tempY;
    //array3_double tempZ;
    //array3_double temp2X;
    //array3_double temp2Y;
    //array3_double temp2Z;
    /*! and some for MaxwellImage */
    array3_double imageX;
    array3_double imageY;
    array3_double imageZ;
    array3_double Dx;
    array3_double Dy;
    array3_double Dz;
    array3_double vectX;
    array3_double vectY;
    array3_double vectZ;
    array3_double divC;
    array3_double arr;

    // *********** SOURCES *************

    // implicit moments obtained from Moments integrator
    //MImoments *miMoments;
    //array3_double rhoc;
    /*! Implicit charge density, defined on central points of the cell */
    //array3_double rhoh;

    // external magnetic field defined on nodes
    //
    array3_double Bx_tot;
    array3_double By_tot;
    array3_double Bz_tot;
    array3_double Bx_ext;
    array3_double By_ext;
    array3_double Bz_ext;

    array3_double Bx_smooth;
    array3_double By_smooth;
    array3_double Bz_smooth;

    // external current, defined on nodes
    // (only used for reporting net current)
    array3_double  *Jx_ext;
    array3_double  *Jy_ext;
    array3_double  *Jz_ext;
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

    /*! CG tolerance criterium for stopping iterations */
    double CGtol;
    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;

    // OpenBC implementation

    injInfoFields *injFieldsLeft, *injFieldsRight, *injFieldsTop, *injFieldsBottom, *injFieldsFront, *injFieldsRear;

    void BoundaryConditionsB(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
      int nx, int ny, int nz);
    void BoundaryConditionsE(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
      int nx, int ny, int nz);
    void BoundaryConditionsEImage(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
      int nx, int ny, int nz);
};

// called by GMRES and CG
void PoissonImage(double *image, double *vector, void** registered_data);

// deprecated
//typedef EMfields3D Field;

#endif
