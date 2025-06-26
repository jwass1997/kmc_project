#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --mem=16gb

echo "$@"

module load compiler/gnu/14.1
module load/boost/1.80.0

"$@"