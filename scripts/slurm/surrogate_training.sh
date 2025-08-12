#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:A100:1
#SBATCH --time=3:00:00
#SBATCH --mem=32G
module load devel/cuda

python3 scripts/python_scripts/surrogate_model.py