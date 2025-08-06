#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:A100:1
#SBATCH --time=04:00:00
#SBATCH --mem=32G

module load devel/cuda
#conda activate py311

python3 scripts/python_scripts/surrogate_model.py