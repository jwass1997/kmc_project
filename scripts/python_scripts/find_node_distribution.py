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

def suggest_spd(trial, 
                prefix, 
                D,
                log_std_bounds=(-3, 3),
                offdiag_bounds=(-1, 1)):
    
    L = np.zeros(shape=(D,D))

    """
    Diagonal elements are constrained (>0)
    """
    for i in range(D):
        log_L_ii = trial.suggest_float(f"{prefix}_logL_{i}{i}", *log_std_bounds)

        L[i, i] = np.exp(log_L_ii)

        """
        Off-diagonal elements are not constrained
        """
        for j in range(i):

            L[i, j] = trial.suggest_float(f"{prefix}_L_{i}{j}", *offdiag_bounds)
    
    return L @ L.T

def suggest_mean_on_disc(trial,
                         prefix,
                         radius):
    
    theta = trial.suggest_float(f"{prefix}_theta", 0.0, 2*np.pi)
    u = trial.suggest_float(f"{prefix}_u", 0.0, 1.0)

    mean_rho = radius*np.sqrt(u)

    mean_x = mean_rho*np.cos(theta)
    mean_y = mean_rho*np.sin(theta)

    return np.array([mean_x, mean_y])

def sample(means,
           covs,
           weights,
           radius,
           N):
    
    positions = []
    cum_sum_weights = np.cumsum(weights)

    while len(positions) < N:

        u = np.random.uniform(0.0, 1.0) * cum_sum_weights[-1]
        component = np.searchsorted(cum_sum_weights, u)

        x = np.random.multivariate_normal(mean=means[component], cov=covs[component])

        if np.sqrt(x[0]**2 + x[1]**2) < radius:

            positions.append(x) 
    
    return np.vstack(positions)

def objective(trial):

    PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
    ROOT = PYTHON_SCRIPT_DIR.parents[2]

    CONFIG_DIR = ROOT / "configs"
    BINARY = ROOT / "build" / "kmc_project"
    DATA_DIR = ROOT / "data" / "experiment_0"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    SLURM_SCRIPT = ROOT / "scripts" / "slurm" / "helix_optuna.sh"

    control_indices = [1, 2, 3, 4, 6, 7]
    input_index = 5
    output_index = 0
    num_of_points = 50
    eq_steps = 1_000
    sim_steps = 1_000
    radius = 150

    samples_per_trial = 10

    print(f"[Trial {trial.number}] Starting trial...")

    """
    Set control voltages
    """
    V_MIN = -1.5
    V_MAX = 1.5

    control_voltages = [-0.2, 1.2, -0.9, 0.8, -0.5, 0.9]#np.random.uniform(V_MIN, V_MAX, size=(len(control_indices, )))
    params = {
        idx: val
        for idx, val in zip(control_indices, control_voltages)
    }

    control_voltage_args = []

    for idx, val in params.items():
        control_voltage_args.append(f"--c_v={idx}={val}")

    """
    Suggest parameters of the gaussian mixture components
    """
    K = 3
    means = []
    covs = []

    raw_weights = []
    for k in range (K):
        means.append(suggest_mean_on_disc(trial=trial, prefix=f"mu_{k}", radius=radius))
        covs.append(suggest_spd(trial=trial, prefix=f"cov_{k}", D=2))

        raw_weight = trial.suggest_float(f"raw_w_{k}", 0.0, 1.0)
        raw_weights.append(raw_weight)

    raw_weight_sum = sum(raw_weights)
    weights = [w / raw_weight_sum for w in raw_weights]

    """
    Sample multiple configuration
    """
    
    NODE_DIST_DIR = DATA_DIR / f"trial_{trial.number}"
    NODE_DIST_DIR.mkdir(parents=True, exist_ok=True)

    for s in range(samples_per_trial):
        X = sample(means, covs, weights, radius, 200)
        file = open(str(NODE_DIST_DIR) + f"/acceptors_trial={trial.number}_sample={s}")

        with open(str(NODE_DIST_DIR) + f"/acceptors_trial={trial.number}_sample={s}", "w") as f:
            for p in X:
                f.write(f"{p[0]}/t{p[1]}\n")                
    """
    Simulation flags
    """
    flags = [
        f"--trial={trial.number}",
        f"--numOfSamples={samples_per_trial}",
        f"--input_idx={input_index}",
        f"--output_idx={output_index}",
        f"--vMin={V_MIN}",
        f"--vMax={V_MAX}",
        f"--numOfPoints={num_of_points}",
        f"--equilibriumSteps={eq_steps}",
        f"--simulationSteps={sim_steps}",
        f"--cfg={str(CONFIG_DIR) + '/config.txt'}",
        f"--donorCfg={str(CONFIG_DIR) + '/donors.txt'}",
        f"--acceptorCfg={str(NODE_DIST_DIR)}"
        f"--electrodeCfg={str(CONFIG_DIR) + '/electrodes.txt'}",
        f"--saveFolderPath={str(DATA_DIR)}",
    ]

    output_file_name = ROOT / "slurm_out" / f"find_node_dist{trial.number}_%j.out"
    slurm_cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SLURM_SCRIPT),
        str(BINARY),
        "findNodeDistribution"     
    ] + flags + control_voltage_args

    subprocess.run(slurm_cmd, capture_output=True, text=True)

    NPZ_FILE = DATA_DIR / f"trial_{trial.number}" / f"curves={trial.number}.npz"

    wait_for_file_creation(NPZ_FILE, 10)

    data = np.load(file=NPZ_FILE)
    curve = data["outputCurrent"]
    vs = np.linspace(V_MIN, V_MAX, num_of_points)
    target = target_function(vs)

    losses += mse(curve, target)#pearson_correlation_coefficient(curve, target)
    
    total_loss = losses / samples_per_trial
    print('\x1b[6;30;42m' + f'{trial.number}] Finished with score: {total_loss}' + '\x1b[0m')

    return total_loss

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