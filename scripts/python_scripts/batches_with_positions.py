import numpy as np
from jobs import slurm_batch_of_independant_states
from pathlib import Path

if __name__ == "__main__":

    batch_size = 500
    num_batches = 200

    for i in range(num_batches):
        
        lhc_seed = np.random.randint(low=0, high=2**10)
        thread_base_seed = np.random.randint(low=0, high=2**10)

        slurm_batch_of_independant_states(
            batch_size=100,
            min_V=-1.5,
            max_V=1.5,
            n_acceptors=200,
            n_electrodes=8,
            n_donors=3,
            radius=150.0,
            nu0=1.0,
            a=20.0,
            T=77,
            energy_disorder=0.01,
            electrode_width=60.0,
            min_hop_distance=1.5,
            max_hop_distance=60.0,
            fem_res=100000,
            dist_type="uniform",
            eps=0.1,
            input_idx=1,
            output_idx=0,
            eq_steps=10_000,
            sim_steps=10_000_000,
            num_of_tasks=100,
            LHCSeed=lhc_seed,
            threadBaseSeed=thread_base_seed,
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/pos_batches_1e7"),
            file_name=Path(f"batch_{i}"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"),
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )
