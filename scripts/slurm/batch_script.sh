#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --mem=32gb

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

"$@"