import numpy as np
from typing import Union
import subprocess
import matplotlib.pyplot as plt
from pathlib import Path

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

def simulate_iv_curve(name: Union[float, str], 
                      save_dir, 
                      cfg_dir, 
                      acc_cfg,
                      don_cfg,
                      ele_cfg,
                      seed: int, 
                      params: dict) -> None:

    BINARY = ROOT / "build" / "kmc_project"

    SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "helix_optuna.sh"

    input_index = 1
    output_index = 0
    num_of_points = 100
    eq_steps = 10_000
    sim_steps = 10_000
    V_MIN, V_MAX = -1.5, 1.5

    params = params

    flags = [
        f"--file_name={name}",
        f"--input_idx={input_index}",
        f"--output_idx={output_index}",
        f"--vMin={V_MIN}",
        f"--vMax={V_MAX}",
        f"--numOfPoints={num_of_points}",
        f"--equilibriumSteps={eq_steps}",
        f"--simulationSteps={sim_steps}",
        f"--cfg={cfg_dir}",
        f"--acceptorCfg={acc_cfg}",
        f"--donorCfg={don_cfg}",
        f"--electrodeCfg={ele_cfg}",
        f"--saveFolderPath={save_dir}",
        f"--seed={seed}"
    ]

    control_voltage_args = []

    for idx, val in params.items():
        control_voltage_args.append(f"--c_v={idx}={val}")

    output_file_name = ROOT / "slurm_out" / "iv_curve_%j.out"
    slurm_cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SLURM_SCRIPT),
        str(BINARY),
        "findControlVoltages"     
    ] + flags + control_voltage_args 

    subprocess.run(slurm_cmd, capture_output=True, text=True)

    #NPZ_FILE = DATA_DIR / f"data_point_{name}.npz"

if __name__ == "__main__":

    voltage_dict = {
        2: 0.5,
        3: -1.2,
        4: 0.0,
        5: -0.6,
        6: 1.3,
        7: -0.9,
        }
    name = "400_acceptors"

    DATA_DIR = ROOT / "data" / "iv_curves"
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ACC_CFG = ROOT / "configs" / "acceptor_vonMises_beta.txt"
    DON_CFG = ROOT / "configs" / "donors.txt"
    ELE_CFG = ROOT / "configs" / "electrodes.txt"
    CFG = ROOT / "configs" / "config.txt"

    simulate_iv_curve(name=name, 
                      seed=65, 
                      save_dir=str(DATA_DIR), 
                      cfg_dir=str(CFG), 
                      acc_cfg=str(ACC_CFG), 
                      don_cfg=str(DON_CFG), 
                      ele_cfg=str(ELE_CFG), 
                      params=voltage_dict)
    