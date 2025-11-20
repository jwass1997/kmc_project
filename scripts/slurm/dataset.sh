#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00:00
#SBATCH --mem=32GB

PYTHON="$HOME/.conda/envs/py311/bin/python"

# (optional) sanity checks
which "$PYTHON" || true
"$PYTHON" -V || true

$PYTHON scripts/python_scripts/create_dataset.py