/*******************************************************************************************
  Particles3D.h  -  Class for particles of the same species, in a 2D space and 3 component velocity
  -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
 ********************************************************************************************/

#ifndef Part2D_H
#define Part2D_H

#include "Particles3Dcomm.h"
//#include "TimeTasks.h"

/**
 * 
 * Class for particles of the same species, in a 2D space and 3 component velocity
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 * I would combine Particles3D and Particles3Dcomm into a
 * single SpeciesPcls class. I believe that the current
 * code actually assumes that Particles3D does not add
 * any data to Particles3Dcomm, as I discovered when
 * I tried to add data to this class. -eaj
 */
class Particles3D:public Particles3Dcomm {

  public:
    /** constructor */
    //Particles3D();
    Particles3D(int species, const Setting& setting_) :
      Particles3Dcomm(species, setting_)
    {}
    /** destructor */
    ~Particles3D(){}
    /** Initial condition: uniform in space and motionless */
    void uniform_background(int is, const_arr4_double rhocs);
    /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
      <ul>
      <li> dim = 0 --> constant velocity on X direction </li>
      <li> dim = 1 --> constant velocity on Y direction </li>
      </ul>
      */
    void constantVelocity(double vel, int dim);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void maxwellian(const_arr4_double rhocs);
    /** Force Free initialization (JxB=0) for particles */
    void force_free(const_arr4_double rhocs);
    /** Rotate velocities in plane XY of angle theta */
    void RotatePlaneXY(double theta);
    /** mover with a Predictor-Corrector Scheme */
    void mover_PC(const_arr4_double aosEMfield);
    /** array-of-structs version of mover_PC */
    void mover_PC_AoS(const_arr4_double aosEMfield);
    /* vectorized version of previous */
    void mover_PC_AoS_vec(const_arr4_double aosEMfield);
    /* mic particle mover */
    void mover_PC_AoS_vec_intr(const_arr4_double aosEMfield);
    /* this computes garbage */
    void mover_PC_AoS_vec_onesort(const_arr4_double aosEMfield);
    /** vectorized version of mover_PC **/
    void mover_PC_vectorized(const_arr4_double aosEMfield);
    /** relativistic mover with a Predictor-Corrector scheme */
    int mover_relativistic(const_arr4_double aosEMfield);
   private:
    /** repopulate particles in a single cell */
    void populate_cell_with_particles(int i, int j, int k, double q,
      double dx_per_pcl, double dy_per_pcl, double dz_per_pcl);
   public:
    /** repopulate particles in boundary layer */
    void repopulate_particles();
    /*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
    double deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center);

#ifdef BATSRUS
    /*! Initial condition: given a fluid model (BATSRUS) */
    void MaxwellianFromFluid(Collective *col);
    /*! Initiate dist. func. for a single cell form a fluid model (BATSRUS) */
    void MaxwellianFromFluidCell(Collective *col, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, longid* ParticleID);
#endif

};

#endif
