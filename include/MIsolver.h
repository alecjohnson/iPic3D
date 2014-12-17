#ifndef MIsolver_h
#define MIsolver_h
// forward declarations
class MImoments;
class Collective;
class Grid3DCU;
class VCtopology3D;
class EMfields3D;
class Kinetics;
class SpeciesMoms;
class Timing;
class OutputWrapperTXT;
class OutputWrapperFPP;
#include "arraysfwd.h"
//
#include "Setting.h"
#include "asserts.h"
#include "assert.h"

namespace iPic3D
{
  // generic class to solve with the moment-implicit method
  class MIsolver
  {
  private:
    const Setting &setting;
    const Collective* col;
    const VCtopology3D* vct;
    const Grid3DCU* grid;
    SpeciesMoms   *speciesMoms;
    MImoments     *miMoments;
    EMfields3D    *EMf; // implicit field solver
    array4_double *fieldForPcls;
    Kinetics      *kinetics;
    //
    // output
    //
    OutputWrapperTXT *outputWrapperTXT;
    OutputWrapperFPP *outputWrapperFPP;
    //
    // convenience variables
    const int ns;
    const int nxn;
    const int nyn;
    const int nzn;
    const int nxc;
    const int nyc;
    const int nzc;

    double tstart;

  public: // accessors
    OutputWrapperFPP& fetch_outputWrapperFPP(){
      assert(outputWrapperFPP);
      return *outputWrapperFPP;
    }
  public:
    ~MIsolver();
    MIsolver(int argc, const char **argv);
    void run();
  protected:
    // virtual so inheriting application can override
    virtual void set_initial_conditions();
    void initialize_output();
    void initialize();
    void init_kinetics_from_restart();
    void init_fields_from_restart();
    void init_from_restart();
    void initDipole();
    void initGEM();
    void initBATSRUS();
    // these should be moved into inheriting class
    void initOriginalGEM();
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation();
    void initGEMDipoleLikeTailNoPert();
    void initGEMnoPert();
    void initRandomField();
    void initForceFree();
  protected:
    void accumulate_moments();
    void compute_moments();
    void advance_Efield();
    void move_particles();
    void advance_Bfield();
    //
    // output methods
    //
    void WriteRestart(int cycle);
    void WriteConserved(int cycle);
    void WriteVelocityDistribution(int cycle);
    void WriteVirtualSatelliteTraces();
    void WriteFields(int cycle);
    void WriteParticles(int cycle);
    void WriteOutput(int cycle);
    void Finalize();

  private:
    void send_field_to_kinetic_solver();
    void set_fieldForPcls();
    int FirstCycle();
    int FinalCycle();
    bool is_rank0();
    void startTiming();
    void stopTiming();
    virtual bool I_am_kinetic_solver(){return true;}
    virtual bool I_am_field_solver(){return true;}

  protected: // accessors
    EMfields3D& fetch_EMfields(){return *EMf;}
    MImoments& fetch_miMoments(){return *miMoments;}
    const MImoments& get_miMoments()const{return *miMoments;}
    SpeciesMoms& fetch_speciesMoms(){return *speciesMoms;}
    const Kinetics& get_kinetics()const{return *kinetics;}
    const Collective& get_col()const{return setting.col();}
    //const Grid3DCU& get_grid()const{return setting.grid();}
    //const VCtopology3D& get_vct()const{return setting.vct();}
    const array4_double& get_fieldForPcls()const{return *fieldForPcls;}
  };
}

#endif // MIsolver_h
