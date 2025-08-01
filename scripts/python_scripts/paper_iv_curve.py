import numpy as np
import matplotlib.pyplot as plt
import subprocess

from pathlib import Path
from scipy.stats import qmc

def create_batch(batch_size,
                 num_of_points,
                 num_of_curves,
                 input_idx,
                 output_idx,
                 eq_steps,
                 sim_steps,
                 seed,
                 batch_id,
                 folder_name,
                 cfg,
                 acceptor_cfg,
                 donor_cfg,
                 electrode_cfg):
    
    PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
    ROOT = PYTHON_SCRIPT_DIR.parents[1]
    #print(ROOT)
    BINARY = ROOT / "build" / "kmc_project"
    SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "batch_script.sh"

    ACC_DIR = ROOT / f"{acceptor_cfg}"
    DON_DIR = ROOT / f"{donor_cfg}"
    ELE_DIR = ROOT / f"{electrode_cfg}"

    CONFIG_DIR = ROOT / f"{cfg}"

    DATA_DIR = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/{folder_name}")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "sbatch",
        f"--output=slurm_out/{batch_id}.out",
        f"--error=slurm_out/{batch_id}.err",
        str(SLURM_SCRIPT),
        str(BINARY),
        "batchRun",
        f"--batchSize={batch_size}",
        f"--numOfPoints={num_of_points}",
        f"--numOfCurves={num_of_curves}",
        f"--inputIdx={input_idx}",
        f"--outputIdx={output_idx}",
        f"--seed={seed}",         
        f"--equilibriumSteps={eq_steps}",
        f"--simulationSteps={sim_steps}",
        f"--cfg={CONFIG_DIR}",
        f"--acceptorCfg={ACC_DIR}",
        f"--donorCfg={DON_DIR}",
        f"--electrodeCfg={ELE_DIR}",
        f"--savePath={DATA_DIR}",
        f"--batchID={batch_id}"
    ]

    subprocess.run(cmd, capture_output=True, text=True)

if __name__ == "__main__":

    electrode_pairs = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 0]
    ]

    batch_size = 1
    num_of_points = 100
    num_of_curves = 1
    eq_steps = 10_000
    sim_steps = 10_000_000

    cfg = "configs/config.txt"
    electrode_cfg = "configs/electrodes.txt"
    folder_name = "iv_data"
    num_devices = 30
    num_of_curves = 1

    for p in range(len(electrode_pairs)):
        for i in range(num_devices):
            seed = np.random.randint(low=0, high=2**31 - 1)
            #print(f"{electrode_pairs[p][0]} {electrode_pairs[p][1]}")
            input_idx = electrode_pairs[p][0]
            output_idx = electrode_pairs[p][1]
            acceptor_cfg = f"configs/iv_curve_configs/acceptors_{i}.txt"
            donor_cfg = f"configs/iv_curve_configs/donors_{i}.txt"

            create_batch(
                batch_size,
                num_of_points,
                num_of_curves,
                input_idx,
                output_idx,
                eq_steps,
                sim_steps,
                seed,
                f"iv_{input_idx}_{output_idx}_{i}",
                folder_name,
                cfg,
                acceptor_cfg,
                donor_cfg,
                electrode_cfg
            )
