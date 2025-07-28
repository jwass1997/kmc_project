import numpy as np
import os
import subprocess

from pathlib import Path

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

BINARY = ROOT / "build" / "kmc_project"
SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "helix_optuna.sh"

ACC_DIR = ROOT / "configs" / "acceptor_normal.txt"
DON_DIR = ROOT / "configs" / "donors.txt"
ELE_DIR = ROOT / "configs" / "electrodes.txt"

DATA_DIR = ROOT / "data"
SWEEP_SAVE_DIR = DATA_DIR / "test_sweep"
SWEEP_SAVE_DIR.mkdir(parents=True, exist_ok=True)

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

input_idx = 1
output_idx = 0
num_of_points = 100
eq_steps = 10_000
sim_steps = 10_000
V_MIN = -1.5
V_MAX = 1.5

num_of_lines = 10
sig_low = 0.01
sig_high = 0.2
t_low = 77.
t_high = 200.

sigmas = np.linspace(sig_low, sig_high, num_of_lines)
temps = np.linspace(t_low, t_high, num_of_lines)

CONFIG_DIR = ROOT / "configs"
SWEEP_CONFIGS = CONFIG_DIR / "t_sig_sweep_configs"
SWEEP_CONFIGS.mkdir(parents=True, exist_ok=True)

cfg_params = {
        "nAcceptors": 200,
        "nDonors": 3,
        "nElectrodes": 8,
        "radius": 150.0,
        "nu0": 1.0,
        "a": 20.0,
        "T": 77.0,
        "energyDisorder": 0.0,
        "electrodeWidth": 60.0,
        "minHopDistance": 1.5,
        "maxHopDistance": 60.0,
        "noDimension": 1
    }

cfg_idx = 0
for t in temps:

    cfg_params["T"] = t
    #print(cfg_params["T"])
    cfg_file = SWEEP_CONFIGS / f"config_{cfg_idx}.txt"
    with cfg_file.open("w") as f:
        for key, val in cfg_params.items():
            f.write(f"{key} {val}\n")

    run_dir = SWEEP_SAVE_DIR / f"run_{cfg_idx}"
    run_dir.mkdir(parents=True, exist_ok=True)
    for i in range(sigmas.shape[0]):
        run_seed = np.random.randint(low=0, high=2**31 - 1)
        flags = [
            f"--file_name=sweep_{i}",
            f"--paramName=sigma",
            f"--paramValue={sigmas[i]}",
            f"--sampleSize={30}",
            f"--input_idx={input_idx}",
            f"--output_idx={output_idx}",
            f"--vMin={V_MIN}",
            f"--vMax={V_MAX}",
            f"--numOfPoints={num_of_points}",
            f"--equilibriumSteps={eq_steps}",
            f"--simulationSteps={sim_steps}",
            f"--cfg={cfg_file}",
            f"--acceptorCfg={str(ACC_DIR)}",
            f"--donorCfg={str(DON_DIR)}",
            f"--electrodeCfg={str(ELE_DIR)}",
            f"--saveFolderPath={run_dir}",
            f"--seed={run_seed + idx}"
        ]

        output_file_name = ROOT / "slurm_out" / "t_sig_sweep_%j.out"

        slurm_CMD = [
            "sbatch",
            #f"--output={output_file_name}",
            str(SLURM_SCRIPT),
            str(BINARY),
            "2DSweep"
        ] + flags + control_voltage_args    

        subprocess.run(slurm_CMD, capture_output=True, text=True)
    
    cfg_idx += 1


print(CONFIG_DIR)
print(SWEEP_CONFIGS)
print(f"ROOT:\n{ROOT}\n")
print(f"PYTHON_SCRIPT_DIR:\n{PYTHON_SCRIPT_DIR}\n")