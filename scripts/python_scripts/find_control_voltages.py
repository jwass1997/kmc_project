import numpy as np
import os
from pathlib import Path
import subprocess
import time
import optuna
import matplotlib.pyplot as plt

def wait_for_file(file_path, timeout=600, check_interval=5):

    start_time = time.time()
    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            raise TimeoutError(f"File {file_path} was not created within {timeout} seconds.")
        time.sleep(check_interval)

def wait_for_file_creation(file_path, check_interval=10):
    
    print(f"Waiting for file: {file_path}")
    while not os.path.exists(file_path):
        time.sleep(check_interval)
    print(f"File found: {file_path}")

def target_function(x: np.array) -> np.array:

    target = 1 / (1 + np.exp(4*x))

    return target

def pearson_correlation_coefficient(x: np.array, y: np.array) -> float:

    N = x.shape[0]

    x_bar = np.mean(x)
    y_bar = np.mean(y)
    std_x = np.std(x, ddof=1)
    std_y = np.std(y, ddof=1)

    x_c = (x - x_bar) / std_x
    y_c = (y - y_bar) / std_y

    pcc = np.sum(x_c*y_c) / (N - 1)

    return pcc   

def mse(x: np.array, y:np.array) -> float:

    N = x.shape[0]

    x_bar = np.mean(x)
    y_bar = np.mean(y)
    std_x = np.std(x, ddof=1)
    std_y = np.std(y, ddof=1)

    x_c = (x - x_bar) / std_x
    y_c = (y - y_bar) / std_y

    mse = np.mean((y_c - x_c)**2)

    return mse

def objective(trial):

    PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
    ROOT = PYTHON_SCRIPT_DIR.parents[1]

    BINARY = ROOT / "build" / "kmc_project"
    CONFIG_DIR = ROOT / "configs"
    DATA_DIR = ROOT / "data" / "trials"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "helix_optuna.sh"

    print(f"[Trial {trial.number}] Starting trial...")

    V_MIN = -1.5
    V_MAX = 1.5

    control_indices = [1, 2, 3, 4, 6, 7]
    input_index = 5
    output_index = 0
    num_of_points = 50
    eq_steps = 10_000
    sim_steps = 10_000

    params = {
        idx: trial.suggest_float(f"c_{idx}", V_MIN, V_MAX)
        for idx in control_indices
    }

    flags = [
        f"--file_name={trial.number}",
        f"--input_idx={input_index}",
        f"--output_idx={output_index}",
        f"--vMin={V_MIN}",
        f"--vMax={V_MAX}",
        f"--numOfPoints={num_of_points}",
        f"--equilibriumSteps={eq_steps}",
        f"--simulationSteps={sim_steps}",
        f"--configs={str(CONFIG_DIR)}",
        f"--saveFolderPath={str(DATA_DIR)}",
    ]

    control_voltage_args = []

    for idx, val in params.items():
        control_voltage_args.append(f"--c_v={idx}={val}")

    output_file_name = ROOT / "slurm_out" / f"find_voltage_trial{trial.number}_%j.out"
    slurm_cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SLURM_SCRIPT),
        str(BINARY),
        "findControlVoltages"     
    ] + flags + control_voltage_args

    subprocess.run(slurm_cmd, capture_output=True, text=True)

    NPZ_FILE = DATA_DIR / f"data_point_{trial.number}.npz"

    wait_for_file_creation(NPZ_FILE, 10)

    data = np.load(file=NPZ_FILE)
    curve = data["outputCurrent"]
    vs = np.linspace(V_MIN, V_MAX, num_of_points)
    target = target_function(vs)

    score = mse(curve, target)#pearson_correlation_coefficient(curve, target)
    print('\x1b[6;30;42m' + f'{trial.number}] Finished with score: {score}' + '\x1b[0m')

    return score

if __name__ == "__main__":

    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(),
        pruner=optuna.pruners.MedianPruner() 
        )

    study.optimize(
        objective,
        n_trials=10,
        timeout=None,
        gc_after_trial=True
        )