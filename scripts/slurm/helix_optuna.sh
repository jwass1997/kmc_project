#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=3:00:00
#SBATCH --mem=64gb

module load compiler/gnu/14.1
module load lib/boost/1.80.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "$@"

"$@"