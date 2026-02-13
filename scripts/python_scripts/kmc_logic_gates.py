import numpy as np
from jobs import slurm_single_device
from pathlib import Path

input_index_0 = 3
input_index_1 = 5
output_idx = 0

num_voltages = 8

GATE = 'NAND'
alpha = 1.0
max_v = 0.5
gate_inputs = [
    [0.0, 0.0],
    [0.0, max_v],
    [max_v, 0.0],
    [max_v, max_v]
]

control_indices = [1, 2, 4, 6, 7]
control_voltages = [-0.38804012537002563, -0.969171404838562, 0.4153963327407837, 0.8823429346084595, -0.9568776488304138]


data_type='uni_1e7'

# CFG_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/config.txt')
# ACC_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/vMB_gradient_away_output.txt')
# DON_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/uniform_donors_1.txt')
# ELE_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/electrodes.txt')

# CFG_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/config.txt')
# ACC_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/vMB_gradient_towards_output.txt')
# DON_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt')
# ELE_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/electrodes.txt')

# CFG_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt')
# ACC_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt')
# DON_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt')
# ELE_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt')

CFG_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/config.txt')
ACC_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/uniform_acceptors_0.txt')
DON_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/uniform_donors_0.txt')
ELE_PATH = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/electrodes.txt')

if __name__ == '__main__':

    voltages = [0.0] * num_voltages

    for pair in gate_inputs:
        voltages[output_idx] = 0.0
        voltages[input_index_0] = pair[0]
        voltages[input_index_1] = pair[1]
        for k, c_i in enumerate(control_indices):
            voltages[c_i] =  control_voltages[k]
        print(control_voltages)
        slurm_single_device(
            control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            control_volts = voltages,
            output_idx=output_idx,
            eq_steps=100_000,
            sim_steps=1_000_000,
            num_intervals=100,
            seed=np.random.randint(low=1, high=2**31 - 1),
            cfg=CFG_PATH,
            acc_cfg=ACC_PATH,
            don_cfg=DON_PATH,
            ele_cfg=ELE_PATH,
            save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_new/logic_gates_data/{data_type}_gate_data_alpha={alpha}"),
            file_name=Path(f"{data_type}_gate={GATE}_{pair[0]}_{pair[1]}_new_1"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/helix_single.sh"),
            OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )   