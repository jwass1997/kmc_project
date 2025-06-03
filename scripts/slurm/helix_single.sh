#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --mem=8gb

cd $SLURM_SUBMIT_DIR
echo "$SLURM_SUBMIT_DIR"
cd build 

./kmc_project singleRun --equilibriumSteps=100000 --simulationSteps=1000000 --deviceName=111