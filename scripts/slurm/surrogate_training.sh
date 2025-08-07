#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:A40:1
#SBATCH --time=01:30:00
#SBATCH --mem=16G

module load devel/cuda
#conda activate py311

python3 scripts/python_scripts/surrogate_model.py