#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:01:00
#SBATCH --mem=8gb

echo "$@"

"$@"