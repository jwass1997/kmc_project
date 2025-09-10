#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:A100:1
#SBATCH --time=3:00:00
#SBATCH --mem=32G

module load devel/cuda

PYTHON="$HOME/.conda/envs/py311/bin/python"

# (optional) sanity checks
which "$PYTHON" || true
"$PYTHON" -V || true
nvidia-smi || true

$PYTHON scripts/python_scripts/train_sm.py \
--data_dir=/gpfs/bwfor/work/ws/hd_gy283-my_data/sm_batches_1e6_mixed_eps=0.8 \
--num_batches=150 \
--batch_size=128 \
--normalize=1 \
--hd=90 \
--num_layers=10 \
--num_epochs=1000 \
--lr=0.001 \
--save_name=sm_mixed_eps=0.8