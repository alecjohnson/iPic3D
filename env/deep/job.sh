#!/bin/bash
#
# This is an example job script for deep.zam.kfa-juelich.de
#
# submit your job with "qsub job.sh"
#
# if you want an interactive job, you can do like this:
#   qsub -I -X -l nodes=2:ppn=16 -l walltime=1:00:00
# This essentially sets the resource limits that are available
# to you in your interactive session.
#
# The PBS lines below give default arguments that would
# otherwise have to be specified on the qsub command line.
# arguments that appear in the qsub command line override
# these directives.
#
# For more information, type "man qsub".
#
# name of job 
#PBS -N ipic3djob
#PBS -l nodes=1:ppn=16
#PBS -v tpt=1
#PBS -l walltime=00:05:00
# oe: merge stderr into stdout; eo: converse; n: separate (default)
#PBS -j oe
# list of users to whom to send email
#PBS -M Alec.Johnson@wis.kuleuven.be
# when to send mail: a=abort,b=beginning,e=end; default: -m a
#PBS -m abe

module load ipic
cd $PBS_O_WORKDIR

# export LD_LIBRARY_PATH=/usr/local/intel/mkl/lib/intel64/:/usr/local/intel/composer_xe_2013.1.117/compiler/lib/intel64/:$LD_LIBRARY_PATH
# # Scalasca configuration
# #export EPK_FILTER=
# #NEXUS="/homeb/zam/knobi/opt/scalasca/1.4.2/bin/scalasca -analyze -f filter.epik"
export NUM_PROCS=16
export NUM_THREADS="-env OMP_NUM_THREADS=$((32/$PROCS))"
rm -rf data/*
mpiexec -np $NUM_PROCS $NUM_THREADS ./iPic3D src/inputfiles/GEM.inp
# 
