#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=05:00:00
#SBATCH --mem=16gb

cd $SLURM_SUBMIT_DIR
echo "$SLURM_SUBMIT_DIR"
cd build 

./kmc_project findControlVoltages \ 
    --name=test \
    --c_v 1=0.5 --c_v 2=-0.5 --c_v 3=1.5 --c_v 4=1.5 --c_v 6=1.5 --c_v 7=1.5 \
    --input_idx=5 \
    --output_idx=0 \
    --vMin=-1.5 \
    --vMax=1.5 \
    --numOfPoints=100 \
    --equilibriumSteps=10000 \
    --simulationSteps=100000 \
    --saveFolderPath=../trials