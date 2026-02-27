import numpy as np
import os
import time
import subprocess
import pandas as pd

from pathlib import Path
from jobs import slurm_single_IV

if __name__ == '__main__':
    data_type = 'grad_towards_1e7'
    func_type = 'Sine'
    input_idx = 4
    maxl2 = 2.0
    seed = 42
    control_indices = [1, 2, 3, 5, 6, 7]

    top_K_num = 25
    noise_strengths = [1, 2, 5, 10, 15, 20]
    if data_type == 'uni_1e7':

        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )

        for k in range(top_K_num):
            for l in range(len(noise_strengths)):
                slurm_single_IV(
                    numOfPoints=50,
                    inputIdx=input_idx,
                    outputIdx=0,
                    control_indices=control_indices,
                    control_volts=df[k],
                    minVoltage=-1.5,
                    maxVoltage=1.5, 
                    eq_steps=100_000,
                    sim_steps=1_000_000,
                    num_intervals=100,
                    seed=np.random.randint(low=1, high=2**30),
                    cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/config.txt"),
                    acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/jiggled_configs/jiggled_{l}.txt"),
                    don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/uniform_donors_0.txt"),
                    ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/electrodes.txt"),
                    save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}"),
                    file_name=Path(f"topk={k}_jiggled={l}"),
                    BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                    SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                    OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                )
    
    if data_type == 'data_vMB_ring_1e7':

        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )

        for k in range(top_K_num):
            for l in range(len(noise_strengths)):
                slurm_single_IV(
                    numOfPoints=50,
                    inputIdx=input_idx,
                    outputIdx=0,
                    control_indices=control_indices,
                    control_volts=df[k],
                    minVoltage=-1.5,
                    maxVoltage=1.5, 
                    eq_steps=100_000,
                    sim_steps=1_000_000,
                    num_intervals=100,
                    seed=np.random.randint(low=1, high=2**30),
                    cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
                    acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/jiggled_configs/jiggled_{l}.txt"),
                    don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
                    ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
                    save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}"),
                    file_name=Path(f"topk={k}_jiggled={l}"),
                    BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                    SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                    OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                )
    
    if data_type == 'grad_towards_1e7':

        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )

        for k in range(top_K_num):
            for l in range(len(noise_strengths)):
                slurm_single_IV(
                    numOfPoints=50,
                    inputIdx=input_idx,
                    outputIdx=0,
                    control_indices=control_indices,
                    control_volts=df[k],
                    minVoltage=-1.5,
                    maxVoltage=1.5, 
                    eq_steps=100_000,
                    sim_steps=1_000_000,
                    num_intervals=100,
                    seed=np.random.randint(low=1, high=2**30),
                    cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/config.txt"),
                    acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/jiggled_configs/jiggled_{l}.txt"),
                    don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt"),
                    ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/electrodes.txt"),
                    save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}"),
                    file_name=Path(f"topk={k}_jiggled={l}"),
                    BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                    SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                    OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                )

    if data_type == 'grad_away_1e7':

        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
        )

        for k in range(top_K_num):
            for l in range(len(noise_strengths)):
                slurm_single_IV(
                    numOfPoints=50,
                    inputIdx=input_idx,
                    outputIdx=0,
                    control_indices=control_indices,
                    control_volts=df[k],
                    minVoltage=-1.5,
                    maxVoltage=1.5, 
                    eq_steps=100_000,
                    sim_steps=1_000_000,
                    num_intervals=100,
                    seed=np.random.randint(low=1, high=2**30),
                    cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/config.txt"),
                    acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/jiggled_configs/jiggled_{l}.txt"),
                    don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/uniform_donors_1.txt"),
                    ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/electrodes.txt"),
                    save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}"),
                    file_name=Path(f"topk={k}_jiggled={l}"),
                    BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                    SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                    OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                )
