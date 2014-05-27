#include "IDgenerator.h"
#include "MPIdata.h"

static int get_rank(){return MPIdata::get_rank();}

void doubleIDgenerator::init(double lowest, double highest)
{
   double num_procs = MPIdata::get_nprocs();
   // technically should add 1 to the difference.
   // in integer case could overflow.
   // after division or floor it doesn't matter anyway.
   double max_pcls = highest - lowest;
   // allocate range of IDs unique to this processor
   double max_pcls_per_proc = floor(max_pcls / num_procs);
   double counter_offset_for_this_proc
     = lowest + max_pcls_per_proc * MPIdata::get_rank();
   // allocate range of IDs unique to each thread
   double max_pcls_per_thread = floor(max_pcls_per_proc / num_threads_in_this_proc);
   for(int i=0;i<num_threads_in_this_proc;i++)
   {
     counter[i] = counter_offset_for_this_proc + i*max_pcls_per_thread;
   }
}

doubleIDgenerator::doubleIDgenerator(
  double lowest, double highest, int num_threads)
{
  if(num_threads <= 0)
    num_threads_in_this_proc = omp_get_max_threads();
  else
    num_threads_in_this_proc = num_threads;

  counter = new double[num_threads_in_this_proc];
  init(lowest,highest);
}

