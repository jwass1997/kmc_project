#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2:00:00
#SBATCH --mem=64gb

set -x

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

seed=$1
row_idx=$2

shift 2

file="row_${row_idx}"
job_seed="$(( seed + row_idx ))"

"$@" \
    --file_name="$file" \
    --seed="$job_seed" \