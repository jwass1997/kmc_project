import numpy as np
import os
import subprocess
import time
import torch
from pathlib import Path
from jobs import slurm_single_device
from bayes_opt import BayesianOptimization

output_index = 0

input_index_0 = 3
input_index_1 = 5

num_indices = 8

control_indices = [
    k for k in range(num_indices) 
    if k not in (input_index_0, input_index_1, output_index)
]

pbounds = {
    f'x_{i}': (-1.5, 1.5) for i in control_indices
}

BINARY_PATH = Path('/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project')
SH_SCRIPT_PATH = Path('/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/helix_single.sh')
OUT_DIR_PATH = Path('/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out')

SAVE_FOLDER_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/BO_vMB_ring_XOR')

CFG_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt')
ACC_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt')
DON_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt')
ELE_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt')

class GateObjective:
    def __init__(self):
        self.iterations = 0

    def __call__(self, **kwargs):

        self.iterations += 1
        
        ITER_PATH = SAVE_FOLDER_PATH / f'iter_{self.iterations}'
        ITER_PATH.mkdir(exist_ok=True, parents=True)
        
        print(f'Iteration: {self.iterations}')
        

        cv_list = [kwargs[f'x_{i}'] for i in control_indices]

        input_pairs = [
            [0.0, 0.0],
            [0.0, .5],
            [.5, 0.0],
            [.5, .5]
        ]

        truth = [
            0,
            1,
            1,
            0
        ]

        file_map = {}
        for l, pair in enumerate(input_pairs):
            voltages = np.zeros(shape=num_indices)
    
            voltages[input_index_0] = pair[0]
            voltages[input_index_1] = pair[1]
    
            for m, c_index in enumerate(control_indices):
                voltages[c_index] = cv_list[m]
    
            full_voltage_list = voltages.tolist()
    
            #print(full_voltage_list)
            FILE_NAME = f'gate_a={pair[0]}_b={pair[1]}_iter={self.iterations}'
            out_file = ITER_PATH / f'{FILE_NAME}.npz'
            file_map[(pair[0], pair[1])] = out_file
    
            slurm_single_device(
                control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
                control_volts=full_voltage_list,
                output_idx = output_index,
                eq_steps=10_000,
                sim_steps=1_000_000,
                num_intervals=100,
                seed=np.random.randint(0, 2 ** 31 - 1),
                cfg=CFG_PATH,
                acc_cfg=ACC_PATH,
                don_cfg=DON_PATH,
                ele_cfg=ELE_PATH,
                save_folder=ITER_PATH,
                file_name=Path(FILE_NAME),
                BINARY=BINARY_PATH,
                SH_SCRIPT=SH_SCRIPT_PATH,
                OUT_DIR=OUT_DIR_PATH
            )
            """args = [
                BINARY_PATH,
                f"singleRun",        
                f"--outputIdx={output_index}",
                f"--eqSteps={10_000}",
                f"--simSteps={100_000}",
                f"--numIntervals={100}",
                f"--seed={31}",
                f"--cfg={str(CFG_PATH)}",
                f"--accCfg={str(ACC_PATH)}",
                f"--donCfg={str(DON_PATH)}",
                f"--eleCfg={str(ELE_PATH)}",
                f"--saveFolder={str(ITER_PATH)}",
                f"--fileName={str(FILE_NAME)}"
            ]
            
            voltage_args = [
                f"--c_v={idx}={volt}"
                for idx, volt in enumerate(full_voltage_list)
            ]

            cmd_args = args + voltage_args
            subprocess.run(args=cmd_args)"""
            
        timeout = 600     
        t0 = time.time()

        while True:
            if all(path.exists() for path in file_map.values()):
                break
            if time.time() - t0 > timeout:
                print('Timeout waiting for npz files')
                return -1e9
            time.sleep(0.5)

        currents = []
        desired = []
        
        for pair, y in zip(input_pairs, truth):
            data = np.load(file_map[(pair[0], pair[1])])
            currents.append(np.array(data['current']).item() * 1e21)
            desired.append(y)
        """  highs = np.array([c for c, t in zip(currents, desired) if t == 1])
        lows = np.array([c for c, t in zip(currents, desired) if t == 0])
        
        margin = highs.min() - lows.max()
        Ion = highs.mean()
        Ioff = lows.mean()
        ratio = Ion / (Ioff + 1e-12)
        score = np.maximum(margin, 0.0) * ratio

        return score """
    
        I_min = min(currents)
        I_max = max(currents)

        delta = I_max - I_min
        
        if abs(delta) < 1e-15:
            print("All currents identical, no contrast")
            return -1e9  # terrible fitness

        norm_currents = [(I - I_min) / delta for I in currents]
        errors = [abs(g - y) for g, y in zip(desired, norm_currents)]
        mean_error = sum(errors) / len(errors)

        F0 = 1.0 - mean_error

        denom = max(abs(I_max), abs(I_min), 1e-15)
        bonus = 0.05 * (delta / denom)
        
        fitness = F0 + bonus
        
        return fitness

"""def run_simulation(**kwargs):

    cv_list = [kwargs[f'x_{i}'] for i in control_indices]

    input_pairs = [
        [0.0, 0.5],
        [0.0, 0.0],
        [0.5, 0.5],
        [0.5, 0.0]
    ]

    for l, pair in enumerate(input_pairs):
        voltages = np.zeros(shape=num_indices)

        voltages[input_index_0] = pair[0]
        voltages[input_index_1] = pair[1]

        for m, c_index in enumerate(control_indices):
            voltages[c_index] = cv_list[m]

        full_voltage_list = voltages.tolist()

        #print(full_voltage_list)
        FILE_NAME = Path(f'gate_')

        slurm_single_device(
            control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            control_volts=full_voltage_list,
            minVoltage=-1.5,
            maxVoltage=1.5,
            eq_steps=10_000,
            sim_steps=1_000_000,
            num_intervals=100,
            seed=42,
            cfg=CFG_PATH,
            acc_cfg=ACC_PATH,
            don_cfg=DON_PATH,
            ele_cfg=ELE_PATH,
            save_folder=SAVE_FOLDER_PATH,
            file_name=FILE_NAME,
            BINARY_PATH=BINARY_PATH,
            SH_SCRIPT=SH_SCRIPT_PATH,
            OUT_DIR=OUT_DIR_PATH
        )"""

if __name__ == '__main__':

    gate_objective = GateObjective()

    optimizer = BayesianOptimization(
        f=gate_objective,
        pbounds=pbounds,
        verbose=2,
        random_state=np.random.randint(0, 2 ** 31 - 1),
    )
    optimizer.maximize(
        init_points=5,
        n_iter=50,
    )

    opt_param_list = []
    target_list = []
    opt_param_list.append(np.array(list(optimizer.max['params'].values())).tolist())
    target_list.append(float(optimizer.max['target']))

    np.save(str(SAVE_FOLDER_PATH / 'results.npy'), {'params': np.stack(opt_param_list), 'targets': np.stack(target_list)})