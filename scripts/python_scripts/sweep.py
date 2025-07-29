import numpy as np
import os
import subprocess

from pathlib import Path

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

BINARY = ROOT / "build" / "kmc_project"
SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "line_sweep.sh"

ACC_DIR = ROOT / "configs" / "acceptor_normal.txt"
DON_DIR = ROOT / "configs" / "donors.txt"
ELE_DIR = ROOT / "configs" / "electrodes.txt"

CONFIG_DIR = ROOT / "configs" / "config.txt"

DATA_DIR = ROOT / "data"
SWEEP_SAVE_DIR = DATA_DIR / "200_sweep"
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

n = 200
const_param_name = "temp"
const_param_min = 77.0
const_param_max = 200.0
const_param = np.linspace(const_param_min, const_param_max, n)

var_param_name = "sigma"
var_param_min = 0.01
var_param_max = 0.2

for i in range(n):
    #ROW_DIR = SWEEP_SAVE_DIR / f"row_{i}"
    #ROW_DIR.mkdir(parents=True, exist_ok=True)
    seed_base = np.random.randint(low=0, high=2**31 - 1)
    cmd = [
        "sbatch",
        f"--job-name=row_{i}",
        str(SLURM_SCRIPT),
        str(i),
        str(seed_base),
        str(BINARY),
        "lineSweep",
        f"--constParamName={const_param_name}",
        f"--varParamName={var_param_name}",
        f"--constParam={const_param[i]}",
        f"--varParamMin={var_param_min}",
        f"--varParamMax={var_param_max}",
        f"--N={n}",
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
        f"--saveFolderPath={str(SWEEP_SAVE_DIR)}",
        *control_voltage_args
    ]

    subprocess.run(cmd, capture_output=True, text=True)