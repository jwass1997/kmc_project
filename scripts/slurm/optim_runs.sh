#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mem=16gb

PYTHON="$HOME/.conda/envs/py311/bin/python"

$PYTHON scripts/python_scripts/run_multiple_optims.py --num_runs 200 --num_iters 100 --optim_type BO --target relu --file_id uni_run_0