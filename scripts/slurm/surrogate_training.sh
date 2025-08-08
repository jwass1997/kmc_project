#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --gres=gpu:A40:1
#SBATCH --time=1:00:00
#SBATCH --mem=16G
module load devel/cuda
export OMP_NUM_THREADS=${SLURM_NTASKS}

python3 scripts/python_scripts/surrogate_model.py