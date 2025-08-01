#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mem=64gb

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

"$@"