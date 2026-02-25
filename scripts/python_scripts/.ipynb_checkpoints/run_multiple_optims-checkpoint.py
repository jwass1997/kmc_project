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

model_dir = f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{args.model}/sm_{args.model}_baseline_0_epochs=1000_dp=0.1_bn=True_hd=90_ls=5_noise=True.pth'
state_dict = torch.load(model_dir, weights_only=False, map_location='cpu')

model = CondSM(state_dict['layer_dims'], state_dict['dropout'], state_dict['batch_norm'])
model.load_state_dict(state_dict['model_state_dict'])
model = model.to(device)
model.eval()

def normalize_over_range(x):
    
    z = (x - x.mean()) / (x.std() + 1e-8) 
    
    return z

def forward_whole_range(model, param_list, input_idx, v_min, v_max, N):

    in_range = torch.linspace(v_min, v_max, N).unsqueeze(1)

    c_volt = torch.tensor(param_list, dtype=torch.float32).expand(N, -1)
    c_volt_l = c_volt[:, :input_idx]
    c_volt_r = c_volt[:, input_idx:]

    input_tensor = torch.cat([c_volt_l, in_range, c_volt_r], dim=1)

    model.eval()
    with torch.no_grad():
        y_out = model(input_tensor).detach().cpu().numpy()

    return y_out

f_target = args.func_type

if f_target == 'ReLU':
    def target_function(x):
    
        m = nn.ReLU()
            
        return m(x)

elif f_target == 'Tanh':
    def target_function(x):
    
        m = nn.Tanh()
            
        return m(x)

elif f_target == 'Sigmoid':
    def target_function(x):
    
        m = nn.Sigmoid()
            
        return m(4*x)

elif f_target == 'Parabola':
    def target_function(x):
            
        return x ** 2

elif f_target == 'Sine':
    def target_function(x):
            
        return torch.sin(3*x)

elif f_target == 'Cubic':
    def target_function(x):
            
        return x ** 3

elif f_target == 'Cosine':
    def target_function(x):
            
        return torch.cos(3*x)

elif f_target == 'Linear':
    def target_function(x):
            
        return x

elif f_target == 'Gaussian':
    def target_function(x):
            
        return torch.exp(-x**2 / 2.0)

elif f_target == 'rbf_gaussian':
    def target_function(x):
        
        return -torch.exp(-1*x**2)

else:
    print('Function to optimize for not found')

def mse_loss_normalized(x, y):

    x = x.float()
    y = y.float()
    x_norm = (x - x.mean()) / (x.std(unbiased=False) + 1e-8)
    y_norm = (y - y.mean()) / (y.std(unbiased=False) + 1e-8)
    
    mse = torch.mean((y_norm - x_norm)**2).item()

    return mse

v_min = -1.5
v_max = 1.5

c_min = -args.c_range
c_max = args.c_range

# --- Number of control electrodes (Normally input dim of SM minus 1) --- #
d = model.layer_dims[0] - 1
pbounds = {
    f'x{i}': (c_min, c_max) for i in range(d) 
}

N = 100
input_range = torch.linspace(v_min, v_max, N, device=device).unsqueeze(1)
input_idx = args.input_index
lam_l1 = args.lam_l1 if args.lam_l1 is not None else 0.0
lam_l2 = args.lam_l2 if args.lam_l2 is not None else 0.0

target = target_function(input_range).to(device)

def f(**kwargs):

    cv_list = [kwargs[f'x{i}'] for i in range(d)]
    c_volt_vec = np.array(cv_list)
    cv_t = torch.tensor(cv_list, dtype=torch.float32, device=device).expand(N, -1)
    cv_t_l = cv_t[:, :input_idx]
    cv_t_r = cv_t[:, input_idx:]

    input_tensor = torch.cat([cv_t_l, input_range, cv_t_r], dim=1)
    model.eval()
    with torch.no_grad():
        y_out = model(input_tensor).detach()

    reg = 0.0
    if args.reg_l1:
        reg += lam_l1 * float(np.sum(np.abs(c_volt_vec))) / (d * np.abs(c_max) + 1e-12)
    if args.reg_l2:
        reg += lam_l2 * float(np.linalg.norm(c_volt_vec, ord=2) / (np.sqrt(d) * (abs(c_max) + 1e-12)))

    loss = mse_loss_normalized(y_out, target) + reg
    
    return -loss

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
        opt_param_vals.append(np.array(list(optimizer.max['params'].values())))
        opt_targets.append(optimizer.max['target'].item())  

    Ts = np.vstack(opt_param_vals)
    df = pd.DataFrame(Ts, columns=[f'x{i}' for i in range(Ts.shape[1])])
    df['target'] = opt_targets
    
    df.to_csv(f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO/bayes_opt_data_type={args.model}_input_idx={args.input_index}_func_type={args.func_type}_lambda_1={args.lam_l1}__lambda_2={args.lam_l2}_c_range={args.c_range}.csv', index=False)