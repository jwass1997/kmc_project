import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from bayes_opt import BayesianOptimization
from models import CondSM

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

parser = argparse.ArgumentParser()

parser.add_argument('--model', type=str)
parser.add_argument('--num_runs', type=int)
parser.add_argument('--num_iters', type=int)
parser.add_argument('--optim_type', type=str)
parser.add_argument('--target', type=str)
parser.add_argument('--file_id', type=str)
args = parser.parse_args()

model_dir = Path(args.model)
state_dict = torch.load(model_dir, weights_only=False, map_location='cpu')

model = CondSM(state_dict['layer_dims'], state_dict['dropout'], state_dict['batch_norm'])
model.load_state_dict(state_dict['model_state_dict'])

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

f_target = args.target

if f_target == 'relu':
    def target_function(x):
    
        m = nn.ReLU()
            
        return m(x)

elif f_target == 'tanh':
    def target_function(x):
    
        m = nn.Tanh()
            
        return m(x)

elif f_target == 'sigmoid':
    def target_function(x):
    
        m = nn.Sigmoid()
            
        return m(x)

elif f_target == 'parabola':
    def target_function(x):
            
        return x ** 2

elif f_target == 'cubic':
    def target_function(x):
            
        return x ** 3

elif f_target == 'rbf_gaussian':
    def target_function(x):
        
        return -np.exp(-1*x**2)

else:
    print('Function to optimize for not found')

def mse_loss(x, y):

    x = x.float()
    y = y.float()
    x_norm = (x - x.mean()) / (x.std(unbiased=False) + 1e-8)
    y_norm = (y - y.mean()) / (y.std(unbiased=False) + 1e-8)
    
    mse = torch.mean((y_norm - x_norm)**2).item()

    return mse

v_min = -1.5
v_max = 1.5

N = 100
input_range = torch.linspace(v_min, v_max, N).unsqueeze(1)

# --- Number of control electrodes (Normally input dim of SM minus 1) --- #
d = model.layer_dims[0] - 1
pbounds = {
    f'x{i}': (v_min, v_max) for i in range(d) 
}

input_idx = 0
target = target_function(input_range)

def f(**kwargs):

    cv_list = list(kwargs.values())
    cv_t = torch.tensor(cv_list, dtype=torch.float32).expand(N, -1)
    cv_t_l = cv_t[:, :input_idx]
    cv_t_r = cv_t[:, input_idx:]

    input_tensor = torch.cat([cv_t_l, input_range, cv_t_r], dim=1)

    model.eval()
    with torch.no_grad():
        y_out = model(input_tensor)

    I_mu = y_out[..., 0].unsqueeze(-1)

    loss = mse_loss(I_mu, target)

    return -loss

def multiple_runs(num_runs, num_iter=10):
    opt_param_list = []
    target_list = []
    for i in range(num_runs):
        optimizer = BayesianOptimization(
            f=f,
            pbounds=pbounds,
            verbose=0,
            random_state=np.random.randint(0, 2**31 - 1),
        )
        optimizer.maximize(
            init_points=10,
            n_iter=num_iter,
        )
        opt_param_list.append(np.array(list(optimizer.max['params'].values())).tolist())
        target_list.append(float(optimizer.max['target']))
        
        print(f'{i+1}/{num_runs} | Params: {opt_param_list[i]}')
    
    opt_param_arr = np.stack(opt_param_list)
    target_arr = np.stack(target_list)

    return opt_param_arr, target_arr

optim_type = args.optim_type

if __name__ == '__main__':
    if optim_type == 'BO':
        params, targets = multiple_runs(num_runs=args.num_runs, num_iter=args.num_iters)
        np.save(f'/gpfs/bwfor/work/ws/hd_gy283-my_data/{optim_type}_{f_target}_{args.file_id}.npy', {'params': params, 'targets': targets})
    else:
        print('Optimization algorithm not found')