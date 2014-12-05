#ifndef MIsolver_h
#define MIsolver_h

class OutputWrapperFPP;

namespace iPic3D
{
  // generic class to solve with the moment-implicit method
  class MIsolver
  {
  private:
    const Setting setting;
    Imoments     *iMoments;
    EMfields3D   *EMf; // implicit field solver
    Kinetics     *kinetics;
    FieldForPcls *fieldForPcls;
    //array4_double fieldForPcls;
    Timing        *my_clock;
    //
    // output
    //
    OutputWrapperTXT *outputWrapperTXT;
    OutputWrapperFPP *outputWrapperFPP;
    //
    // convenience variables
    const int ns;

    // accessors

  public:
    EMfields3D& fetch_fieldSolver() {return *EMf;}
    Kinetics& fetch_kinetics() {return *kinetics;}
    
    //
    OutputWrapperFPP& fetch_outputWrapperFPP(){
      assert(outputWrapperFPP);
      return *outputWrapperFPP;
    }
  public:
    ~MIsolver();
    MIsolver(int argc, const char **argv);
    // virtual so inheriting application can override
    virtual void set_initial_conditions();
    int Init();
    void computeMoments();
    void advanceEfield();
    void moveParticles();
    void advanceBfield();
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

    int FirstCycle() { return setting.col().get_first_cycle(); }
    int FinalCycle() { return setting.col().get_final_cycle(); }
    bool is_rank0() { return setting.vct().is_rank0(); }

  public: // accessors
    const Collective& get_col()const{return setting.col()}
    const Grid& get_grid()const{return setting.grid();};
    const VirtualTopology3D& get_vct()const{return setting.vct();}

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
