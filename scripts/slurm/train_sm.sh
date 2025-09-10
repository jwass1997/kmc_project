#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:A100:1
#SBATCH --time=36:00:00
#SBATCH --mem=128G

module load devel/cuda

PYTHON="$HOME/.conda/envs/py311/bin/python"

# (optional) sanity checks
which "$PYTHON" || true
"$PYTHON" -V || true
nvidia-smi || true

$PYTHON scripts/python_scripts/train_sm.py