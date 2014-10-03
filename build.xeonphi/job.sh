#!/bin/bash
#
# template SLURM script for use with MPI
# rush this script with: sbatch job.sh
# monitor jobs with:
#   showq
#    squeue -u `whoami`
# kill jobs with:
#   scancel <jobId>
#
# alternatively you can submit an interactive
# job on a single node with
#
#   srun --pty -n 16 -t 01:00:00 -p development /bin/bash -l
#
# or equivalently: idev -pe 16 1 -m 60 # 1 = single node
#
# job name
#SBATCH -J ipic3d
##SBATCH -o ipic3d.%j.out
##SBATCH -e ipic3d.%j.err
# queue (development, normal, or normal-mic)
#SBATCH -p normal-mic
# number of nodes, not cores (16 cores per node)
#TemplateLine: #SBATCH -N $NUM_NODES
#SBATCH -N 4
# total number of MPI tasks (if omitted, n=N)
#TemplateLine: #SBATCH -n $NUM_PROCS
#SBATCH -n 4
# maximum time
#SBATCH -t 00:05:00
#SBATCH --mail-user=Alec.Johnson@wis.kuleuven.be
#SBATCH --mail-type=ALL  

#module load ipic
#module list

#TemplateLine: export XLEN=$XLEN
export XLEN=2
#TemplateLine: export YLEN=$YLEN
export YLEN=2
# export OMP_NUM_THREADS=8
# NUM_PROCS = nodes * ppn / OMP_NUM_THREADS
#TemplateLine: NUM_PROCS=$NUM_PROCS
NUM_PROCS=$(($XLEN*$YLEN))
DATA=data
#export PATH="/home1/02456/alec/bin:/home1/02456/alec/scripts:$PATH"

#export IPIC_HOME=/home1/02456/alec/ipic3d
#export IPIC_ENV=stampede
#. $IPIC_HOME/env/bashrc

#module use ~/ipic3d/env
#module load ipic-common
#module use ~/ipic3d/env/stampede
#module list
#module load ipic
#module list
set -x
pwd
# generate lists of machines for offloading
scontrol show hostname | tee spawnfile | sed 's/$/-mic0/' > machinefile
#export I_MPI_MIC=1
#export DAPL_UCM_REP_TIME=8000
#export DAPL_UCM_RTU_TIME=4000
#export DAPL_UCM_RETRY=10
#export I_MPI_HYDRA_BRANCH_COUNT=6000
#export I_MPI_OFA_ADAPTER_NAME=mlx4_0
#export I_MPI_DAPL_PROVIDER_LIST="ofa-v2-mlx4_0-1u,ofa-v2-scif0,ofa-v2-mcm-1"
#mpiexec.hydra -iface br0 -genv I_MPI_HYDRA_BRANCH_COUNT 6000 -genv I_MPI_FALLBACK 0 -genv I_MPI_FABRICS shm:dapl -genv OMP_NUM_THREADS 16 -np 4 -machinefile machinefile -env LD_LIBRARY_PATH $MIC_LD_LIBRARY_PATH exec/iPic3D ../inputfiles/GEM.inp
#
# set environment for ibrun.symm -m
#export MIC_PPN=1
#export MIC_OMP_NUM_THREADS=4
#export MIC_MY_NSLOTS=4
# use ibrun for MPI codes, not mpirun or srun
mpiexec.hydra -iface br0 -np $NUM_PROCS -machinefile machinefile -env LD_LIBRARY_PATH $MIC_LD_LIBRARY_PATH exec/iPic3D ../inputfiles/GEM.inp
#ibrun.symm -m "exec/iPic3D ../inputfiles/GEM.inp" -iface br0 -np 4 -machinefile machinefile -env LD_LIBRARY_PATH $MIC_LD_LIBRARY_PATH
#ibrun.symm -m "exec/iPic3D ../inputfiles/GEM.inp" -np 4 -machinefile machinefile
#ibrun.symm -m -np 4 exec/iPic3D ../inputfiles/GEM.inp | tee out.${XLEN}x${YLEN}.txt
