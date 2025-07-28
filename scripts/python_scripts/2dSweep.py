import numpy as np
import os
import subprocess

from pathlib import Path

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

BINARY = ROOT / "build" / "kmc_project"
SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "test.sh"

ACC_DIR = ROOT / "configs" / "acceptor_normal.txt"
DON_DIR = ROOT / "configs" / "donors.txt"
ELE_DIR = ROOT / "configs" / "electrodes.txt"

CONFIG_DIR = ROOT / "configs" / "config.txt"

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
num_of_points = 50
eq_steps = 10_000
sim_steps = 10_000
V_MIN = -1.5
V_MAX = 1.5

n = 100

for i in range(n):
    ROW_DIR = SWEEP_SAVE_DIR / f"row_{i}"
    ROW_DIR.mkdir(parents=True, exist_ok=True)
    seed_base = np.random.randint(low=0, high=2**31 - 1)
    cmd = [
        "sbatch",
        f"--job-name=sweep_row_{i}",
        f"--array=0-{n-1}%100",
        str(SLURM_SCRIPT),
        str(i),
        str(n),
        str(seed_base),
        str(BINARY),
        "2DSweep",
        f"--sampleSize={5}",
        f"--input_idx={input_idx}",
        f"--output_idx={output_idx}",
        f"--vMin={V_MIN}",
        f"--vMax={V_MAX}",
        f"--numOfPoints={num_of_points}",
        f"--equilibriumSteps={eq_steps}",
        f"--simulationSteps={sim_steps}",
        f"--cfg={str(CONFIG_DIR)}",
        f"--acceptorCfg={str(ACC_DIR)}",
        f"--donorCfg={str(DON_DIR)}",
        f"--electrodeCfg={str(ELE_DIR)}",
        f"--saveFolderPath={str(ROW_DIR)}",
        *control_voltage_args
    ]

    subprocess.run(cmd, capture_output=True, text=True)