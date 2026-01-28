import numpy as np
from pathlib import Path
import time
from jobs import slurm_batch_from_single_state

batch_size = 2000

CFG = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/config.txt")
ACC_CFG = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/uniform_acceptors_0.txt")
DON_CFG = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/uniform_donors_0.txt")
ELE_CFG = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/electrodes.txt")

if __name__ == '__main__':
    for i in range(100):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=100_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**30 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=CFG,
            acc_cfg=ACC_CFG,
            don_cfg=DON_CFG,
            ele_cfg=ELE_CFG,
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7"),
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )
    
    """ for i in range(100):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/vMB_gradient_towards_output.txt"),
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt"),
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7"),
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """