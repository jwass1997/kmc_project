#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --job-name=test_gpu
#SBATCH --output=logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:H200:1
#SBATCH --time=6:00:00
#SBATCH --mem=32gb

module load devel/cuda

PYTHON="$HOME/.conda/envs/py311/bin/python"

srun $PYTHON /home/hd/hd_hd/hd_gy283/kmc_project/scripts/python_scripts/flow_training.py
