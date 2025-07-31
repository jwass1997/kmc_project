#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=5:00:00
#SBATCH --mem=64gb
#SBATCH --output=slurm_out/pixel_%A_%a.out
#SBATCH --error=slurm_out/pixel_%A_%a.err

set -x

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

seed=$1
const_param_min=$2
const_param_max=$3
N=$4

shift 4

row_idx=${SLURM_ARRAY_TASK_ID}

const_param_step=$(echo "scale=6; ($const_param_max - $const_param_min)/($N - 1)" | bc -l)

const_param=$(echo "scale=6; $const_param_min + $row_idx * $const_param_step" | bc -l)

file="row_${row_idx}"
job_seed="$(( seed + row_idx ))"

"$@" \
    --file_name="$file" \
    --seed="$job_seed" \
    --constParam="$const_param"