#ifndef _Parameters_h_
#define _Parameters_h_

// use of namespace provides a more flexible, succinct
// singleton via "using Parameters" but in this case might
// be better called GlobalParams.  Parameters seems to
// be the appropriate name for the "Collective" class.
//
namespace Parameters
{
  enum Enum
  {
    null=0,
    SoA, // struct of arrays
    AoS, // array of structs
    // for moments type
    AoSvec,
    SoAvec,
    AoSintr,
    // for mover type
    SoA_vec_onesort,
    AoS_vec_onesort,
    SoA_vec_resort,
    AoS_vec_resort,
    //
    // velocity capping
    //
    // velocity cap scale
    fraction_of_subdomain,
    fraction_of_domain,
    absolute_velocity,
    fraction_of_light_speed,
    // velocity cap method
    box_capped,
    sphere_capped,
    //
    // pusher
    //
    relativistic,
    classical,
  };

  void init_parameters();

  bool get_USING_AOS();
  bool get_SORTING_SOA();
  bool get_SORTING_PARTICLES();
  // for resorting particles with each iteration of mover
  bool get_RESORTING_PARTICLES();
  inline bool get_USING_XAVG() { return get_RESORTING_PARTICLES(); }
  bool get_VECTORIZE_MOMENTS();
  //bool get_VECTORIZE_MOVER();
  Enum get_MOVER_TYPE();
  Enum get_MOMENTS_TYPE();

  int get_multiple_of_vector_width_in_doubles();
  // blocksize and numblocks for use in BlockCommunicator
  int get_blockSize();
  int get_numBlocks();
  bool get_doWriteOutput();
  //
  int get_num_smoothings();
  bool use_perfect_smoothing();
  bool use_correct_smoothing();
  inline bool use_original_smoothing(){return !use_correct_smoothing();}

  // pusher
  //
  int vel_cap_method();
  int vel_cap_scale_ref();
  double vel_cap_scaled_value();
  int pusher();
}
#endif
