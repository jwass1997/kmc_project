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
        f"--output=slurm_out/batch.out",
        f"--error=slurm_out/batch.err",
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

    batch_size = 1
    num_of_points = 100
    num_of_curves = 30
    input_idx = 1
    output_idx = 0
    eq_steps = 10_000
    sim_steps = 100_000
    batch_id = 80
    seed = np.random.randint(low=0, high=2**30 - 1)
    folder_name = "iv_pair=01"

    acceptor_cfg = "configs/acceptors.txt"
    donor_cfg = "configs/donors.txt"
    electrode_cfg = "configs/electrodes.txt"
    cfg = "configs/config.txt"
        
    create_batch(
        batch_size,
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
        electrode_cfg
    )