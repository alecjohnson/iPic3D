#!/bin/sh
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/apps/intel/13/composer_xe_2013_sp1.1.106/compiler/lib/intel64

# Set number of OpenMP threads per spawned MPI proc 
export OMP_NUM_THREADS=4

../build.xeon/exec/iPic3D_fields $1
