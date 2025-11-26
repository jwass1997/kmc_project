#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mem=16gb

PYTHON="$HOME/.conda/envs/py311/bin/python"

$PYTHON scripts/python_scripts/run_multiple_optims.py --num_runs 20 --num_iters 200 --optim_type BO --target cubic --file_id vMB_0