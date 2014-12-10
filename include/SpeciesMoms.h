#ifndef SpeciesMoms_h
#define SpeciesMoms_h
#include "Alloc.h"

// Moments that kinetcs module writes to
// and fields solver can be ignorant of.
//
class Setting;
class Particles3Dcomm;

// class to accumulate node-centered species moments
// for a single species.
// 
class Moments10
{
  private:
    arr4_double arr;
    int nx;
    int ny;
    int nz;
  public:
    void set_to_zero();

    // fetch accessors (write access)
    arr4_double fetch_arr() { return arr; }

    Moments10(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      arr (nxn, nyn, nzn,10)
    {
    };
    ~Moments10(){};
};

// The ten primitive particle moments of each species
// needed by the implicit moment method
//
// (This code was previously embedded in the EMfields class.)
//
class SpeciesMoms
{
    const Setting& setting;
    // convenience copies
    const int nxc;
    const int nyc;
    const int nzc;
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

  public:
    // accessors

    arr4_double fetch_rhons(){return rhons;}
    const_arr4_double get_rhons()const{return rhons;}

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

    /*** accessors that require computing ***/

    // get current for species in all cells except ghost
    //
    arr3_double ret_Jxsc(int is);
    arr3_double ret_Jysc(int is);
    arr3_double ret_Jzsc(int is);

    ~SpeciesMoms();
    SpeciesMoms(const Setting& setting_);

  // --- operations to compute species moments
  private:

    void setZeroSpeciesMoms(int is);

    /*! sum moments (interp_P2G) versions */
    void sumMoments(const Particles3Dcomm& pcls);
    void sumMoments_vec(const Particles3Dcomm& pcls);
    void sumMoments_AoS(const Particles3Dcomm& pcls);
    void sumMoments_AoS_intr(const Particles3Dcomm& pcls);
    void sumMoments_vectorized(const Particles3Dcomm& pcls);
    void sumMoments_vectorized_AoS(const Particles3Dcomm& pcls);
    void sumMomentsOld(const Particles3Dcomm& pcls);

    /*! adjust densities on boundaries that are not periodic */
    void adjustNonPeriodicDensities(int is);

  public:
    void accumulateMoments(const int is, Particles3Dcomm& pcls);
    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns);

  // operations to calculate implicit moments
  public:
    void calculateJhat(
      arr3_double Jxh,
      arr3_double Jyh,
      arr3_double Jzh,
      const_arr3_double Bx,
      const_arr3_double By,
      const_arr3_double Bz)const;
};

//arr3_double SpeciesMoms::ret_rhocs(int is) { return ca.get_N2C_no_ghosts(is, rhons); }
//arr3_double SpeciesMoms::ret_Jxsc(int is) { return ca.get_N2C_no_ghosts(is, Jxs); }
//arr3_double SpeciesMoms::ret_Jysc(int is) { return ca.get_N2C_no_ghosts(is, Jys); }
//arr3_double SpeciesMoms::ret_Jzsc(int is) { return ca.get_N2C_no_ghosts(is, Jzs); }

#endif // SpeciesMoms_h
