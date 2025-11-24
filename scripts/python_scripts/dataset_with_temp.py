import numpy as np
import time
from pathlib import Path
from jobs import slurm_batch_temp_dependent

num_batches = 200
batch_size = 1_000

if __name__ == '__main__':
    for i in range(200):
        time.sleep(0.1)
        slurm_batch_temp_dependent(
            batch_size=1000,
            min_V=-1.5,
            max_V=1.5,
            Tmax=300.0,
            Tmin=70.0,
            n_acceptors=200,
            n_electrodes=8,
            n_donors=3,
            radius=150.0,
            nu0=1.0,
            a=20.0,
            energy_disorder=0.01,
            electrode_width=60.0,
            min_hop_distance=1.5,
            max_hop_distance=60.0,
            Nr=257,
            Nt=1440,
            dist_type="uniform",
            n_comps=8,
            output_idx=0,
            eq_steps=10_000,
            sim_steps=1_000_000,
            num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31),
            threadBaseSeed=np.random.randint(low=0, high=2**31),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/dataset_wTemp"),
            file_name=f"batch_{i}",
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"),
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )
