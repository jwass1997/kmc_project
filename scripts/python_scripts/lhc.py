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
                 folder_name):
    
    PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
    ROOT = PYTHON_SCRIPT_DIR.parents[1]

    BINARY = ROOT / "build" / "kmc_project"
    SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "batch_script.sh"

    ACC_DIR = ROOT / "configs" / "acceptors.txt"
    DON_DIR = ROOT / "configs" / "donors.txt"
    ELE_DIR = ROOT / "configs" / "electrodes.txt"

    CONFIG_DIR = ROOT / "configs" / "config.txt"

    DATA_DIR = ROOT / "data" / f"{folder_name}"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    voltage_dict = {
    2: 0.5,
    3: -1.2,
    4: 0.0,
    5: -0.6,
    6: 1.3,
    7: -0.9
    }

    control_voltage_args = []

    for idx, val in voltage_dict.items():
        control_voltage_args.append(f"--c_v={idx}={val}")
    
    cmd = [
        "sbatch",
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

    """ sampler = qmc.LatinHypercube(d=7)
    sample = sampler.random(n=5)
    l_bounds = [-1.5] * 7
    u_bounds = [1.5] * 7
    scaled_sample = qmc.scale(sample, l_bounds, u_bounds)
    print(scaled_sample) """

    batch_size = 10
    num_of_points = 100
    num_of_curves = 5
    input_idx = 1
    output_idx = 0
    eq_steps = 10_000
    sim_steps = 100_000
    batch_id = 112
    seed = np.random.randint(low=0, high=2**31 - 1)
    folder_name = "test_batch"
        
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
        folder_name
    )