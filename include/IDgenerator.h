#ifndef IDgenerator_h
#define IDgenerator_h

#include "ompdefs.h"

// Class to generate unique double-precision IDs
//
// Increasing the number of unique IDs to 2^64 would be possible
// using uint64_t or via large double precision IDs created by
// grouping processors e.g. into 2048 different groups each
// having a unique exponent.
//
// Idea of how to assign unique double precision IDs
// using most of the double precision bits:
//
//   MPI_Comm_rank returns int, so assume that
//   there are at most 2^31 processors.
//   Putting 2^30 processors in 2^10 groups
//   means that there are 2^20 processors in a group.
//   let group_idx range from 1 to 2^10.
//   let ID_increment = group_idx.
//   set counter as if there are only 2^20 processors.
//   then multiply counter by group_idx.
//
// But I need to be careful about the sign bits
// in the mantissa and exponent.
//
// We can work out the details when we really need more
// than 2^53 unique particle IDs.
//
class doubleIDgenerator
{
  double * counter;
  int num_threads_in_this_proc;
  // largest consecutive positive integer representable by double
  static const double DINTMAX = 0x20000000000000p0; // 2^53
 private:
  void init(double lowest, double highest);
 public:
  doubleIDgenerator(double lowest=0, double highest=DINTMAX, int num_threads=-1);
  ~doubleIDgenerator() { delete [] counter; }
 public: 
  double getID()
  {
    return counter[omp_get_thread_num()]++;
  }
};

#endif // IDgenerator_h
