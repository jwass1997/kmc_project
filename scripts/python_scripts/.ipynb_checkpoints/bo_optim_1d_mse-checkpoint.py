import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
from pathlib import Path
from bayes_opt import BayesianOptimization
from models import CondSM

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str)
parser.add_argument('--num_init_points', type=int)
parser.add_argument('--num_runs', type=int)
parser.add_argument('--num_iters', type=int)
parser.add_argument('--input_index', type=int)
parser.add_argument('--func_type', type=str)
parser.add_argument('--reg_l1', action='store_true', help='enable L1 regularization')
parser.add_argument('--lam_l1', type=float)
parser.add_argument('--reg_l2', action='store_true', help='enable L2 regularization')
parser.add_argument('--lam_l2', type=float)
parser.add_argument('--c_range', type=float)
args = parser.parse_args()

model_dir = (
    f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{args.model}/'
    f'sm_{args.model}_baseline_0_epochs=1000_dp=0.1_bn=True_hd=90_ls=5_noise=True.pth'
)
state_dict = torch.load(model_dir, weights_only=False, map_location='cpu')

model = CondSM(state_dict['layer_dims'], state_dict['dropout'], state_dict['batch_norm'])
model.load_state_dict(state_dict['model_state_dict'])
model = model.to(device).eval()

# ----------------------------
# Target function
# ----------------------------
f_target = args.func_type

if f_target == 'ReLU':
    def target_function(x): return nn.ReLU()(x)
elif f_target == 'Tanh':
    def target_function(x): return nn.Tanh()(x)
elif f_target == 'Sigmoid':
    def target_function(x): return nn.Sigmoid()(4 * x)
elif f_target == 'Parabola':
    def target_function(x): return x ** 2
elif f_target == 'Sine':
    def target_function(x): return torch.sin(3 * x)
elif f_target == 'Cubic':
    def target_function(x): return x ** 3
elif f_target == 'rbf_gaussian':
    def target_function(x): return -torch.exp(-1 * x**2)
else:
    raise ValueError('Function to optimize for not found')

# ----------------------------
# Affine-aligned loss
# ----------------------------
def affine_aligned_mse(y, t, eps=1e-12, enforce_positive_gain=False):
    """
    Fits a,b to minimize ||a*y + b - t||^2 and returns MSE residual.

    y: [N] or [N,1] or [B,N]
    t: [N] or [N,1] or [B,N]
    returns: mse (scalar if single curve), plus a,b (scalars)
    """
    # squeeze possible trailing dim
    if y.ndim == 2 and y.shape[1] == 1:
        y = y.squeeze(1)
    if t.ndim == 2 and t.shape[1] == 1:
        t = t.squeeze(1)

    if y.ndim == 1:
        y = y.unsqueeze(0)  # [1,N]
    if t.ndim == 1:
        t = t.unsqueeze(0).expand_as(y)

    y_mean = y.mean(dim=1, keepdim=True)
    t_mean = t.mean(dim=1, keepdim=True)

    y0 = y - y_mean
    t0 = t - t_mean

    var_y = (y0 * y0).mean(dim=1, keepdim=True)
    cov_yt = (y0 * t0).mean(dim=1, keepdim=True)

    a = cov_yt / (var_y + eps)
    if enforce_positive_gain:
        a = torch.clamp(a, min=0.0)

    b = t_mean - a * y_mean

    y_fit = a * y + b
    mse = ((y_fit - t) ** 2).mean(dim=1)  # [B]
    return mse.squeeze(0), a.squeeze(0).squeeze(0), b.squeeze(0).squeeze(0)

# ----------------------------
# Ranges / bounds
# ----------------------------
v_min, v_max = -1.5, 1.5
c_min = -args.c_range
c_max = args.c_range

# d = number of controls = input_dim - 1
d = model.layer_dims[0] - 1
pbounds = {f'x{i}': (c_min, c_max) for i in range(d)}

N = 100
input_range = torch.linspace(v_min, v_max, N, device=device).unsqueeze(1)  # [N,1]
input_idx = args.input_index

lam_l1 = args.lam_l1 if args.lam_l1 is not None else 0.0
lam_l2 = args.lam_l2 if args.lam_l2 is not None else 0.0

target = target_function(input_range).to(device)  # [N,1] typically

def f(**kwargs):
    cv_list = [kwargs[f'x{i}'] for i in range(d)]
    c_volt_vec = np.array(cv_list, dtype=float)

    # build input tensor [N, d+1]
    cv_t = torch.tensor(cv_list, dtype=torch.float32, device=device).expand(N, -1)  # [N,d]
    cv_t_l = cv_t[:, :input_idx]
    cv_t_r = cv_t[:, input_idx:]
    input_tensor = torch.cat([cv_t_l, input_range, cv_t_r], dim=1)  # [N, d+1]

    with torch.no_grad():
        y_out = model(input_tensor)  # [N] or [N,1]

    # affine-aligned loss (best gain+offset per evaluation)
    loss_aff, a_star, b_star = affine_aligned_mse(
        y_out, target,
        enforce_positive_gain=True  # set False if you want to allow negative gain
    )
    loss = float(loss_aff.item())

    # regularization on controls (unchanged)
    reg = 0.0
    if args.reg_l1:
        reg += lam_l1 * float(np.sum(np.abs(c_volt_vec))) / (d * np.abs(c_max) + 1e-12)
    if args.reg_l2:
        reg += lam_l2 * float(np.linalg.norm(c_volt_vec, ord=2) / (np.sqrt(d) * (abs(c_max) + 1e-12)))

    loss_total = loss + reg

    # maximize target => return negative loss
    return -loss_total

if __name__ == '__main__':
    num_iters = args.num_iters
    num_init_points = args.num_init_points
    num_runs = args.num_runs

    opt_param_vals = []
    opt_targets = []

    for i in range(num_runs):
        print("Run iteration:", i)
        optimizer = BayesianOptimization(
            f=f,
            pbounds=pbounds,
            verbose=1,
            random_state=np.random.randint(0, 2**31 - 1),
        )
        optimizer.maximize(
            init_points=num_init_points,
            n_iter=num_iters,
        )
        print('run done:', i)

        opt_param_vals.append(np.array(list(optimizer.max['params'].values()), dtype=float))
        opt_targets.append(float(optimizer.max['target']))

    Ts = np.vstack(opt_param_vals)
    df = pd.DataFrame(Ts, columns=[f'x{i}' for i in range(Ts.shape[1])])
    df['target'] = opt_targets

    out_path = (
        f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse/'
        f'bayes_opt_affine_data_type={args.model}_input_idx={args.input_index}_func_type={args.func_type}_'
        f'lambda_1={args.lam_l1}__lambda_2={args.lam_l2}_c_range={args.c_range}.csv'
    )
    df.to_csv(out_path, index=False)
    print("Saved:", out_path)