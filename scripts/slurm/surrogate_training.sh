#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:A100:1
#SBATCH --time=3:00:00
#SBATCH --mem=32G

module load devel/cuda

# Sanity checks (shows which Python is used and if CUDA is visible)
echo "Python: $(which python)"
python - <<'PY'
import sys, torch
print("sys.executable:", sys.executable)
print("torch:", torch.__version__)
print("CUDA available?:", torch.cuda.is_available())
print("Torch CUDA version:", getattr(torch.version, "cuda", None))
PY

conda activate py311

python3 scripts/python_scripts/surrogate_model.py \
--data_dir=/gpfs/bwfor/work/ws/hd_gy283-my_data/sm_batches_1e6_mixed_eps=0.2 \
--num_batches=150 \
--batch_size=128 \
--normalize=1 \
--hd=90 \
--num_layers=10 \
--num_epochs=1000 \
--lr=0.001 \
--save_name=sm_mixed_eps=0.2