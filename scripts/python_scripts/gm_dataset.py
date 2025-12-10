import numpy as np
from pathlib import Path
import time
from jobs import slurm_batch_from_single_state

batch_size = 2000

if __name__ == '__main__':
    """for i in range(200):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**30 - 1), threadBaseSeed=np.random.randint(low=0, high=2**30 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/gm_data/configs/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/gm_data/configs/gm_acceptors_0.txt"),
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/gm_data/configs/uniform_donors_0.txt"),
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/gm_data/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/gm_data"),
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )"""
    
    for i in range(100):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/uni_data/configs/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/uni_data/configs/uniform_acceptors_0.txt"),
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/uni_data/configs/uniform_donors_0.txt"),
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/uni_data/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/paper_datasets/uni_data"),
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )