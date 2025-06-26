#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=05:00:00
#SBATCH --mem=16gb

module load compiler/gnu/14.1
module load lib/boost/1.80.0

echo "$@"

"$@"