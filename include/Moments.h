#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

class Particles3Dcomm;

// moments that must be communicated from the kinetic solver
// to the field solver of the moment-implicit method.
class Imoments
{
  private: // parameters
    // access to global data
    const Setting& setting;
    // dimensions of arrays
    const int nx, ny, nz;
    // number of species
    const int ns;
  private: // moments
    // densities of each species
    vector_array3_double rhons;

    // implicit current density defined on nodes
    //
    array3_double Jxh;
    array3_double Jyh;
    array3_double Jzh;

  public:
    const vector_array3_double& get_rhons()const{return rhons;}
    const_arr3_double get_Jxh()const{return Jxh;}
    const_arr3_double get_Jyh()const{return Jyh;}
    const_arr3_double get_Jzh()const{return Jzh;}

  public:
    Imoments(const Setting& setting_)
    : setting(setting_),
      nxn(setting.grid().getNXN()),
      nyn(setting.grid().getNYN()),
      nzn(setting.grid().getNZN()),
      ns(setting.col().getNs()),
      rhons(ns, nxn, nyn, nzn),
      Jxh(nxn, nyn, nzn),
      Jyh(nxn, nyn, nzn),
      Jzh(nxn, nyn, nzn)
    {
    }
};

// class to accumulate node-centered species moments
// for a single species
// 
class Moments10
{
  private:

    //arr4_double arr;
    int nx;
    int ny;
    int nz;
  public:
    void set_to_zero();

    // fetch accessors (write access)
    //arr4_double fetch_arr() { return arr; }

    Moments10(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      //arr (nxn, nyn, nzn,10)
    {
    };
    ~Moments10(){};
};

// The ten primitive particle moments of each species
// needed by the implicit moment method
//
// (This code was previously embedded in the EMfields class.)
//
class Pmoments
{
    const Setting& setting;
    const int nxn;
    const int nyn;
    const int nzn;
    const int ns;

    // *** species-specific moments needed to construct implicit moments ***

    /*! charge density on central points of cells */
    array4_double rhons;
    //
    // species-specific current densities defined on nodes
    //
    array4_double Jxs;
    array4_double Jys;
    array4_double Jzs;
    //
    // pressure tensor components, defined on nodes
    //
    array4_double pXXsn;
    array4_double pXYsn;
    array4_double pXZsn;
    array4_double pYYsn;
    array4_double pYZsn;
    array4_double pZZsn;

    // *** temporary arrays for calculations ***

    /* temporary arrays for summing moments */
    int sizeMomentsArray;
    Moments10 **moments10Array;

    /*! fetch array for summing moments of thread i */
    Moments10& fetch_moments10Array(int i){
      assert_le(0,i);
      assert_lt(i,sizeMomentsArray);
      return *(moments10Array[i]);
    }
    int get_sizeMomentsArray() { return sizeMomentsArray; }

    // *** moments for output ***

    // These should be computed by the class responsible
    // for output when output is performed, not computed
    // in this class in case they might be needed.
    //
    /*! Charge density, defined on central points of the cell */
    //array3_double rhoc;
    /*! Charge density, defined on nodes */
    //array3_double rhon;
    //
    // net current density defined on nodes
    //
    //array3_double Jx;
    //array3_double Jy;
    //array3_double Jz;

  public:
    // accessors

    arr4_double fetch_rhons(){return rhons;}
    const_arr4_double get_rhons()const{return rhons;}
    //arr4_double fetch_rhocs(){return rhocs;}

    // deprecated
    //
    //double getRHOns(int X,int Y,int Z,int is)const{return rhons.get(is,X,Y,Z);}
    //double*** getRHOns(int is){return rhons.fetch_arr4()[is];}
    //arr4_double getRHOns(){return rhons;}
    //
    arr4_double getJxs() { return Jxs; }
    arr4_double getJys() { return Jys; }
    arr4_double getJzs() { return Jzs; }
    double getJxs(int X,int Y,int Z,int is)const{return Jxs.get(is,X,Y,Z);}
    double getJys(int X,int Y,int Z,int is)const{return Jys.get(is,X,Y,Z);}
    double getJzs(int X,int Y,int Z,int is)const{return Jzs.get(is,X,Y,Z);}
    arr4_double getpXXsn() { return pXXsn; }
    arr4_double getpXYsn() { return pXYsn; }
    arr4_double getpXZsn() { return pXZsn; }
    arr4_double getpYYsn() { return pYYsn; }
    arr4_double getpYZsn() { return pYZsn; }
    arr4_double getpZZsn() { return pZZsn; }

    // accessors for moments for reporting

    // densities per species:
    //
    //double getRHOcs(int X,int Y,int Z,int is)const{return rhocs.get(is,X,Y,Z);}
    /* density on cells without ghost cells */
    //arr3_double getRHOcs(int is);

    //arr3_double getRHOc() { return rhoc; }
    //arr3_double getRHOn() { return rhon; }
    //double getRHOc(int X, int Y, int Z) const { return rhoc.get(X,Y,Z);}
    //double getRHOn(int X, int Y, int Z) const { return rhon.get(X,Y,Z);}

    /*** accessors that require computing ***/

    // get current for species in all cells except ghost
    //
    arr3_double ret_Jxsc(int is);
    arr3_double ret_Jysc(int is);
    arr3_double ret_Jzsc(int is);

  Pmoments(const Setting& setting_)
    setting(setting_),
    nxn(setting_.grid().get_nxn()),
    nyn(setting_.grid().get_nyn()),
    nzn(setting_.grid().get_nzn()),
    ns(setting_.col().getNs()),
    //
    // species-specific quantities
    //
    rhons (ns, nxn, nyn, nzn),
    Jxs   (ns, nxn, nyn, nzn),
    Jys   (ns, nxn, nyn, nzn),
    Jzs   (ns, nxn, nyn, nzn),
    pXXsn (ns, nxn, nyn, nzn),
    pXYsn (ns, nxn, nyn, nzn),
    pXZsn (ns, nxn, nyn, nzn),
    pYYsn (ns, nxn, nyn, nzn),
    pYZsn (ns, nxn, nyn, nzn),
    pZZsn (ns, nxn, nyn, nzn),

    // temporary array for reporting cell-centered values
    arr (nxc-2,nyc-2,nzc-2),

    // only used for initialization and output
    //rhocs (ns, nxc, nyc, nzc),
  {}

  private:

    // operations to calculate implicit moments

    /*! Calculate the three components of Pi(implicit pressure)
        cross image vector */
    void PIdot(arr3_double PIdotX, arr3_double PIdotY, arr3_double PIdotZ,
      const_arr3_double vectX, const_arr3_double vectY, const_arr3_double vectZ, int ns);

    /*! Calculate rho hat, Jx hat, Jy hat, Jz hat */
    void calculateHatFunctions();

    // operations to calculate moments

    /*! set to 0 primary moments */
    void setZeroPrimaryMoments();
    /*! set to 0 all densities derived from primary moments */
    void setZeroDerivedMoments();
    /*! Sum rhon over species */
    void sumOverSpecies();
    /*! Sum current over different species */
    void sumOverSpeciesJ();

    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns);
    //
    /*! sum moments (interp_P2G) versions */
    void sumMoments(const Particles3Dcomm* part);
    void sumMoments_vec(const Particles3Dcomm* part);
    void sumMoments_AoS(const Particles3Dcomm* part);
    void sumMoments_AoS_intr(const Particles3Dcomm* part);
    void sumMoments_vectorized(const Particles3Dcomm* part);
    void sumMoments_vectorized_AoS(const Particles3Dcomm* part);
    void sumMomentsOld(const Particles3Dcomm& pcls);

    /*! adjust densities on boundaries that are not periodic */
    void adjustNonPeriodicDensities(int is);

};

#endif
