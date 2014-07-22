#!/bin/bash
#
# template SLURM script for use with MPI
# rush this script with: sbatch job.sh
# monitor jobs with:
#   showq
#    squeue -u `whoami`
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
#SBATCH -o ipic3d.%j.out
##SBATCH -e ipic3d.%j.err
# queue
#SBATCH -p development
# number of nodes, not cores (16 cores per node)
#SBATCH -N 2
# total number of MPI tasks (if omitted, n=N)
#SBATCH -n 32
# maximum time
#SBATCH -t 00:05:00
##SBATCH --mail-user=Alec.Johnson@wis.kuleuven.be
##SBATCH --mail-type=ALL  

module load ipic
module list

export XLEN=4
export YLEN=4
# export OMP_NUM_THREADS=8
# NUM_PROCS = nodes * ppn / OMP_NUM_THREADS
NUM_PROCS=$(($XLEN*$YLEN))
echo "running on $NUM_PROCS cpus ..."

# use ibrun for MPI codes, not mpirun or srun
ibrun -np $NUM_PROCS ./iPic3D GEM.inp | tee out.${XLEN}x${YLEN}.txt
