#!/bin/bash
#SBATCH --partition=serial

#SBATCH --qos=84c-1d_serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3G
#SBATCH --job-name="Dizon_test"
#SBATCH --output=Dizon_test.%J.out

#SBATCH --error=Dizon_test.%J.err
#SBATCH --mail-type=ALL
#SBATCH --requeue

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# Set stack size to unlimited
ulimit -s unlimited

# MAIN
./PEOA
