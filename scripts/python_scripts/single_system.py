import os
from pathlib import Path
import time
import subprocess

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

BINARY = ROOT / "build" / "kmc_project"
CONFIG_DIR = ROOT / "configs"
DATA_DIR = ROOT / "data"

flags = [
    f"--cfg={str(CONFIG_DIR) + '/config.txt'}",
    f"--acceptorCfg={str(CONFIG_DIR) + '/acceptor_normal.txt'}",
    f"--donorCfg={str(CONFIG_DIR) + '/donors.txt'}",
    f"--electrodeCfg={str(CONFIG_DIR) + '/electrodes.txt'}",
    f"--save_path={str(DATA_DIR)}",
    f"--equilibriumSteps={10000}",
    f"--simulationSteps={100000}",
    f"--deviceName={999}"
]

SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "helix_single.sh"

if __name__ == "__main__":

    #print(ROOT)
    #print(str(PYTHON_SCRIPT_DIR))
    #print(CONFIG_DIR)
    #print(DATA_DIR)
    #print(SLURM_SCRIPT)

    slurm_cmd = [
        "sbatch", 
        str(SLURM_SCRIPT),
        str(BINARY),
        "singleRun",
    ] + flags
    #print(slurm_cmd)
    subprocess.run(slurm_cmd, capture_output=True, text=True)