#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:15:00
#SBATCH --mem=8gb

echo "$@"

"$@"