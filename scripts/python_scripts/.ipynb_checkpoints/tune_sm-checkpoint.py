import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import optuna

from sm import NeuralNet

print(optuna.__version__)
optuna.logging.set_verbosity(optuna.logging.WARNING)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

def find_control_voltages_1D(model, target, N, input_idx, min_voltage, max_voltage, n_startup_trials=1000, n_trials=2000, PRINT_EVERY=100):
    
    model.eval()

    def mse_loss(x: torch.Tensor, y: torch.Tensor):

        N = x.shape[0]
        x = x.float()
        y = y.float()
        x_norm = (x - x.mean()) / (x.std(unbiased=False) + 1e-8)
        y_norm = (y - y.mean()) / (y.std(unbiased=False) + 1e-8)
        
        mse = torch.mean((y_norm - x_norm)**2).item()

        return mse

    v_input = torch.linspace(min_voltage, max_voltage, N).unsqueeze(1)

    def objective(trial):

        cv_list = [trial.suggest_float(f"x_{i}", min_voltage, max_voltage) for i in range(6)]
        cv_tensor = torch.tensor(cv_list, dtype=torch.float32)
        cv_input = cv_tensor.expand(N, -1)
        
        cv_input_left = cv_input[:, :input_idx]
        cv_input_right = cv_input[:, input_idx:]

        input_tensor = torch.cat([cv_input_left, v_input, cv_input_right], dim=-1)
        with torch.no_grad():
            y_pred = model(input_tensor).detach().cpu()
    
        mse = mse_loss(y_pred, target)

        if trial.number % PRINT_EVERY == 0:
            print('\x1b[6;30;42m' + f'{trial.number}] Finished with score: {mse: .6f}' + '\x1b[0m')

        return mse
    
    tpe = optuna.samplers.TPESampler(multivariate=True, group=True, n_startup_trials=n_startup_trials)
    study = optuna.create_study(direction="minimize", sampler=tpe)
    study.optimize(objective, n_trials=n_trials)

    return study

def find_control_voltages_2D(model, target_flat, N, input_indices, min_voltage, max_voltage, n_startup_trials=500, n_trials=1000, PRINT_EVERY=1000):
    
    model.eval()
    device = next(model.parameters()).device
    print(device)
    in_features = 7
    i1, i2 = sorted(input_indices)
    target_flat = target_flat.to(device).float()

    v = torch.linspace(min_voltage, max_voltage, N, device=device)
    v1, v2 = torch.meshgrid(v, v, indexing="ij")
    v1_input = v1.ravel()
    v2_input = v2.ravel()

    cv_cols = [c for c in range(in_features) if c not in (i1, i2)]

    def mse_loss(x, y):
        # flatten to 1D and normalize globally (like the 1D version)
        x = x.reshape(-1).float()
        y = y.reshape(-1).float()
        x = (x - x.mean()) / (x.std(unbiased=False) + 1e-8)
        y = (y - y.mean()) / (y.std(unbiased=False) + 1e-8)
        return torch.mean((y - x) ** 2).item()
    
    def objective(trial):

        cv_list = [trial.suggest_float(f"x_{i}", min_voltage, max_voltage) for i in range(in_features - 2)]
        cv_tensor = torch.tensor(cv_list, dtype=torch.float32, device=device)
        
        input_tensor = torch.empty(N*N, in_features, dtype=torch.float32, device=device)
        for col, val in zip(cv_cols, cv_tensor):
            input_tensor[:, col] = val
        
        input_tensor[:, i1] = v1_input
        input_tensor[:, i2] = v2_input

        with torch.no_grad():
            y_pred = model(input_tensor)

        mse = mse_loss(y_pred, target_flat)      

        if trial.number % PRINT_EVERY == 0:
            print('\x1b[6;30;42m' + f'{trial.number}] Finished with score: {mse: .6f}' + '\x1b[0m')

        return mse
    
    tpe = optuna.samplers.TPESampler(multivariate=True, group=True, n_startup_trials=n_startup_trials)
    study = optuna.create_study(direction="minimize", sampler=tpe)
    study.optimize(objective, n_trials=n_trials)

    return study

if __name__ == "__main__":

    def target_function(t: torch.Tensor):
        N = t.shape[0]
        return -0.7*torch.sin(3*t) + torch.sigmoid(7*t)
    
    min = -1.5
    max = 1.5
    N = 100
    x = torch.linspace(min, max, N)
    target = target_function(x).unsqueeze(1)

    state_dict = torch.load(f="/home/hd/hd_hd/hd_gy283/kmc_project/models/SM_1e6_vonMises_beta_0.pth", map_location=torch.device("cpu"))
    #print(state_dict)

    in_features = 7
    out_features = 1
    hidden_dim = 90
    num_layers = 10
    
    model = NeuralNet(in_features=in_features, out_features=out_features, hidden_dim=hidden_dim, num_layers=num_layers)
    model.load_state_dict(state_dict=state_dict)

    study = find_control_voltages_1D(model=model, target=target, N=N, min_voltage=min, max_voltage=max, input_idx=1)

    print("Best value:", study.best_value)
    print("Best params:", study.best_params)    