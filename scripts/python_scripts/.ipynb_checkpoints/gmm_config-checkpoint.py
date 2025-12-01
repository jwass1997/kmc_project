import numpy as np
import time
from pathlib import Path
from jobs import slurm_batch_from_single_state

num_batches = 200
batch_size = 1_500

if __name__ == '__main__':

    for i in range(num_batches):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=1_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_fewer_electrodes/cfgs/config.txt"), 
            acc_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_fewer_electrodes/cfgs/acceptors.txt"), 
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_fewer_electrodes/cfgs/donors.txt"), 
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_fewer_electrodes/cfgs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )