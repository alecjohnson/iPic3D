
#include <mpi.h>
#include <stdarg.h>
#include "TimeTasks.h"
#include "asserts.h"
#include "MPIdata.h" // for get_rank
#include "parallel.h"
#include "debug.h"
#include "Collective.h"

/** implementation of declarations in utility/TimeTasks.h **/

TimeTasks timeTasks;

static const char *taskNames[] = // order must agree with Tasks in TimeTasks.h
{
  "none",
  "moments",
  "fields",
  "particles",
  "last",
  "bfield",
  "moment_pcl_sorting",
  "moment_accumulation",
  "moment_reduction",
  "mover_pcl_sorting",
  "mover_pcl_moving",
  "transpose_pcls_to_AoS",
  "transpose_pcls_to_SoA",
  "write_fields",
  "write_particles",
  "number_of_tasks"
};

const char* TimeTasks::get_taskname(int arg)
{
  assert_le(arg,NUMBER_OF_TASKS);
  return taskNames[arg];
}

void TimeTasks::resetCycle()
{
  for(int e=0;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]=0.;
    communicate[e]=0.;
    active[e]=false;
    stack_depth[e]=0;
    start_times[e]=0.;
  }
  active_task=NONE;
  communicating=false;
}
void TimeTasks::start_main_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(is_exclusive(taskid));
  assert_ne(active_task, taskid);
  active_task = taskid;
  assert(!active[taskid]);
  active[taskid]=true;
}
void TimeTasks::start_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(!is_exclusive(taskid));
  assert(!active[taskid]);
  active[taskid]=true;
}
// have to manage the task stack explicitly
void TimeTasks::start_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  if(stack_depth[taskid]==0)
  {
    start_times[taskid]=start_time;
    start_task(taskid);
  }
  stack_depth[taskid]++;
}
void TimeTasks::end_main_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  end_task(taskid, start_time);
  active_task = NONE;
}
void TimeTasks::end_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(active[taskid]);
  double now = MPI_Wtime();
  // compute time spent on task
  task_duration[taskid] += now - start_time;
  active[taskid] = false;
}
// have to manage the task stack explicitly
void TimeTasks::end_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  stack_depth[taskid]--;
  assert_ge(stack_depth[taskid],0);
  if(stack_depth[taskid]==0)
  {
    end_task(taskid, start_times[taskid]);
  }
}
void TimeTasks::end_communicating(double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(active_task);
  assert(communicating);
  double additional_communication_time = MPI_Wtime()-start_time;
  communicate[active_task] += additional_communication_time;
  communicating=false;
}

void TimeTasks::print_cycle_times(int cycle,
  double* reported_task_duration,
  double* reported_communicate,
  const char* timing_prefix)
{

  // restrict output to master process
  //
  if(MPIdata::get_rank()) return;
  FILE* file = stdout;
  {
    //fflush(file);
    //fprintf(file,"=== times for cycle %d for rank %d === \n",
    //  cycle,
    //  MPIdata::get_rank());
    fprintf(file, "%stotal  comput commun task\n", timing_prefix);
    for(int e=NONE+1; e<LAST; e++)
    {
      fprintf(file, "%s%6.3f %6.3f %6.3f %s\n",
      timing_prefix,
      reported_task_duration[e], //get_time(e),
      reported_task_duration[e]-reported_communicate[e], //get_compute(e),
      reported_communicate[e], //communicate[e],
      get_taskname(e));
    }

    // report total times
    //
    // get total time spent on exclusive tasks
    //
    double total_task_duration = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total_task_duration += task_duration[i];
    }
    // get total time spent in exclusive tasks spent communicating
    //
    double total_communicate = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total_communicate += communicate[i];
    }
    const double total_computing_time = total_task_duration - total_communicate;
    fprintf(file, "%s%6.3f %6.3f %6.3f %s\n",
      timing_prefix,
      total_task_duration,
      total_computing_time,
      total_communicate,
      "[total times]");

    fprintf(file, "%stime   subtask\n", timing_prefix);
    for(int e=LAST+1; e<NUMBER_OF_TASKS; e++)
    {
      // do not show tasks that are not executed
      double elapsed_time = reported_task_duration[e];
      if(!elapsed_time)
        continue;

      assert_eq(stack_depth[e],0);
      fprintf(file, "%s%6.3f %s\n",
      timing_prefix,
      elapsed_time,
      get_taskname(e));
    }
    
    fflush(file);
  }
}

void TimeTasks::print_cycle_times_min(int cycle)
{
  // assume that only main thread is active
  assert(!omp_get_thread_num());
  // perform all-reduce to get max times for all processes
  double reported_task_duration[NUMBER_OF_TASKS];
  double reported_communicate[NUMBER_OF_TASKS];
  {
    MPI_Allreduce(task_duration,reported_task_duration,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(communicate,reported_communicate,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  }
  print_cycle_times(cycle, reported_task_duration, reported_communicate, "min_| ");
}

void TimeTasks::print_cycle_times_max(int cycle)
{
  // assume that only main thread is active
  assert(!omp_get_thread_num());
  // perform all-reduce to get max times for all processes
  double reported_task_duration[NUMBER_OF_TASKS];
  double reported_communicate[NUMBER_OF_TASKS];
  {
    MPI_Allreduce(task_duration,reported_task_duration,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(communicate,reported_communicate,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }
  print_cycle_times(cycle, reported_task_duration, reported_communicate, "max_| ");
}

void TimeTasks::print_cycle_times_avg(int cycle)
{
  // assume that only main thread is active
  assert(!omp_get_thread_num());
  // perform all-reduce to average times for all processes
  double reported_task_duration[NUMBER_OF_TASKS];
  double reported_communicate[NUMBER_OF_TASKS];
  {
    double total_task_duration[NUMBER_OF_TASKS];
    double total_communicate[NUMBER_OF_TASKS];
    MPI_Allreduce(task_duration,total_task_duration,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(communicate,total_communicate,
      NUMBER_OF_TASKS,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // compute averages.
    const int nprocs = MPIdata::get_nprocs();
    for(int i=0;i<NUMBER_OF_TASKS;i++)
    {
      reported_task_duration[i] = total_task_duration[i]/nprocs;
      reported_communicate[i] = total_communicate[i]/nprocs;
    }
  }
  print_cycle_times(cycle, reported_task_duration, reported_communicate, "avg_| ");
}

void TimeTasks::print_cycle_times(int cycle)
{
  assert(!omp_get_thread_num());

  if(!MPIdata::get_rank())
    fflush(stdout);
  if(!MPIdata::get_rank())
    printf("=== times for cycle %d (main process) ===\n", cycle);
  print_cycle_times(cycle, task_duration, communicate, "main| ");
  if(!MPIdata::get_rank())
    printf("=== times for cycle %d (maximum over all processes) ===\n", cycle);
  print_cycle_times_max(cycle);
  if(!MPIdata::get_rank())
    printf("=== times for cycle %d (averaged over all processes) ===\n", cycle);
  print_cycle_times_avg(cycle);
  if(!MPIdata::get_rank())
    printf("=== times for cycle %d (minimum over all processes) ===\n", cycle);
  print_cycle_times_min(cycle);
}

// The following three methods provide for a hack by which
// the timeTasks copies of all threads are averaged.
// 
void TimeTasks::operator/=(int num)
{
  assert(false); // this method is not in use.
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]/=num;
    start_times[e]/=num;
    communicate[e]/=num;
  }
}
void TimeTasks::operator+=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  communicating = arg.communicating;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    active[e] = arg.active[e];
    task_duration[e]+=arg.task_duration[e];
    stack_depth[e] = arg.stack_depth[e];
    start_times[e]+=arg.start_times[e];
    communicate[e]+=arg.communicate[e];
  }
}
void TimeTasks::operator=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  communicating = arg.communicating;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    active[e] = arg.active[e];
    task_duration[e]=arg.task_duration[e];
    stack_depth[e] = arg.stack_depth[e];
    start_times[e]=arg.start_times[e];
    communicate[e]=arg.communicate[e];
  }
}

TimeTasks_caller_to_set_main_task_for_scope::
TimeTasks_caller_to_set_main_task_for_scope(TimeTasks::Tasks _task) :
  task(_task)
{
  //if(omp_get_thread_num()) return;
  // assume that only one thread is active
  assert(!omp_get_thread_num());
  // Since we only report the time for the main process and
  // main thread, we synchronize when measuring times.  A
  // better solution would be to average the times for all
  // processes and threads.
  //MPI_Barrier(MPI_COMM_WORLD);
  start_time = MPI_Wtime();
  timeTasks.start_main_task(task);
}
TimeTasks_caller_to_set_main_task_for_scope::
~TimeTasks_caller_to_set_main_task_for_scope()
{
  //if(omp_get_thread_num()) return;
  // assume that only one thread is active
  assert(!omp_get_thread_num());

  //MPI_Barrier(MPI_COMM_WORLD);
  timeTasks.end_main_task(task, start_time);
}

TimeTasks_caller_to_set_task_for_scope::
TimeTasks_caller_to_set_task_for_scope(TimeTasks::Tasks _task)
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  task = _task;
  already_active = timeTasks.is_active(task);
  if(!already_active)
  {
    //#pragma omp barrier
    start_time = MPI_Wtime();
    timeTasks.start_task(task);
  }
}
TimeTasks_caller_to_set_task_for_scope::
~TimeTasks_caller_to_set_task_for_scope()
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  if(already_active)
  {
    assert(timeTasks.is_active(task));
  }
  else
  {
    //#pragma omp barrier
    timeTasks.end_task(task, start_time);
  }
}

TimeTasks_caller_to_set_communication_mode_for_scope::
TimeTasks_caller_to_set_communication_mode_for_scope()
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  already_communicating = timeTasks.get_communicating();
  if(!already_communicating)
  {
    start_time = MPI_Wtime();
    timeTasks.set_communicating(true);
  }
}
TimeTasks_caller_to_set_communication_mode_for_scope::
~TimeTasks_caller_to_set_communication_mode_for_scope()
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  if(!already_communicating)
  {
    timeTasks.end_communicating(start_time);
  }
}
