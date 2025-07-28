#!/bin/bash
#SBATCH --output=slurm_out/pixel_%A_%a.out
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --mem=16gb

param1_min=77.0
param1_max=200.0
param1_step=$(echo "scale=6; ($param1_max - $param1_min)/($n - 1)" | bc -l)

param2_min=0.01
param2_max=0.2
param2_step=$(echo "scale=6; ($param2_max - $param2_min)/($n - 1)" | bc -l)

