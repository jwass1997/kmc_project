#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=10:00:00
#SBATCH --mem=64gb

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd $SLURM_SUBMIT_DIR
echo "$SLURM_SUBMIT_DIR"
cd build 

./kmc_project batchRun --batchSize=100 --equilibriumSteps=10000 --simulationSteps=1000000 --batchName=222
