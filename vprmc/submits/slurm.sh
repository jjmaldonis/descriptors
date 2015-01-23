#!/bin/sh

#SBATCH --job-name=vprmc                  # job name
#SBATCH --partition=univ                # default "univ" if not specified
#SBATCH --error=job.%J.err              # error file
#SBATCH --output=job.%J.out             # output file

#SBATCH --time=7-00:00:00               # run time in days-hh:mm:ss

#SBATCH --nodes=1                      # number of nodes requested (n)
#SBATCH --ntasks=1                     # required number of CPUs (n)
#SBATCH --ntasks-per-node=1             # default 16 (Set to 1 for OMP)
#SBATCH --cpus-per-task=1              # default 1 (Set to 16 for OMP)
##SBATCH --mem=16384                    # total RAM in MB, max 64GB  per node
##SBATCH --mem-per-cpu=4000              # RAM in MB (default 4GB, max 8GB)

# set OMP_NUM_THREADS to the number of --cpus-per-task we asked for
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

##SBATCH --export=ALL

echo "Date:"
date '+%s'
echo "Using ACI / HCP / Slurm cluster."
echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Using $SLURM_NODELIST nodes."
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

# Executable
mpirun rmc $SLURM_JOB_ID vp_indexes_t3_top90.txt ../testmodels/Zr50Cu35Al15_10k_start.xyz


echo "Finished on:"
date '+%s'
