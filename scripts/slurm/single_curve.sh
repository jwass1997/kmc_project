#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:5:00
#SBATCH --mem=8gb

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

"$@"