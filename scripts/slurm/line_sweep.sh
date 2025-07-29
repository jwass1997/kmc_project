#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --mem=32gb
#SBATCH --output=slurm_out/pixel_%A_%a.out
#SBATCH --error=slurm_out/pixel_%A_%a.err

set -x

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

row_idx=$1
seed=$2

shift 2

file="row_${row_idx}"
job_seed="$(( seed + row_idx ))"

"$@" \
    --file_name="$file" \
    --seed="$job_seed" 