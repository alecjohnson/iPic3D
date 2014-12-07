#include "Setting.h"
#ifndef MIsolver_h
#define MIsolver_h

class OutputWrapperFPP;
class Imoments;
class EMfields3D;
class Kinetics;

namespace iPic3D
{
  // generic class to solve with the moment-implicit method
  class MIsolver
  {
  private:
    const Setting &setting;
    Imoments      *iMoments;
    EMfields3D    *EMf; // implicit field solver
    Kinetics      *kinetics;
    array4_double *fieldForPcls;
    //array4_double fieldForPcls;
    Timing        *my_clock;
    //
    // output
    //
    OutputWrapperTXT *outputWrapperTXT;
    OutputWrapperFPP *outputWrapperFPP;
    //
    // convenience variables
    //const int ns;

    // accessors

  public:
    //EMfields3D& fetch_fieldSolver() {return *EMf;}
    //Kinetics& fetch_kinetics() {return *kinetics;}
    
    //
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
  protected:
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
    int FirstCycle() { return setting.col().get_first_cycle(); }
    int FinalCycle() { return setting.col().get_final_cycle(); }
    bool is_rank0() { return setting.vct().is_rank0(); }

  protected: // accessors
    EMfields3D& fetch_EMfields(){return *EMf;}
    Imoments& fetch_iMoments(){return *iMoments;}
    Pmoments& fetch_pMoments(){return pMoments;}
    const Kinetics& get_kinetics()const{return *kinetics;}
    const Collective& get_col()const{return setting.col()}
    const Grid& get_grid()const{return setting.grid();};
    const VirtualTopology3D& get_vct()const{return setting.vct();}
  private:
    int Init();
  };

  // solver that chooses initial and boundary
  // conditions based on configuration
  // (move this into apps directory).
  //
  class c_Solver: public MIsolver
  {
    // this class overrides this method
    virtual void set_initial_conditions();
  };
}

#endif MIsolver_h
