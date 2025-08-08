#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --gres=gpu:A40:2
#SBATCH --time=5:00:00
#SBATCH --mem=64gb
module load devel/cuda
export OMP_NUM_THREADS=${SLURM_NTASKS}

python3 scripts/python_scripts/surrogate_model.py