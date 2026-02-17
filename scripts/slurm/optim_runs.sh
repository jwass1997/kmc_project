#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --mem=16gb

PYTHON="$HOME/.conda/envs/py311/bin/python"

$PYTHON scripts/python_scripts/run_multiple_optims.py --model uni_1e7 --num_init_points 5 --num_runs 20 --num_iters 100 --input_idx 3 --func_type Sine --lam 0.0