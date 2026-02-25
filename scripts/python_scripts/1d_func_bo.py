import numpy as np
import os
import time
import subprocess
import pandas as pd

from pathlib import Path
from jobs import slurm_single_IV

if __name__ == "__main__":
  data_type = 'grad_away_1e7'
  func_type = 'Cubic'
  input_idx = 4
  maxl2 = 2.0
  seed = 42
  control_indices = [1, 2, 3, 5, 6, 7]

  df = pd.read_csv(
    f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse/bayes_opt_affine_data_type={data_type}'
    f'_input_idx=3_func_type={func_type}_lambda_1=0.0__lambda_2=1.0_c_range=1.5.csv'
    )

  x_cols = [c for c in df.columns if c.startswith("x")]
  best = df.nsmallest(1, "target")   # or nlargest(1, "target")
  control_volts = best[x_cols].iloc[0].to_list()

  if data_type == 'uni_1e7':
    slurm_single_IV(
      numOfPoints=50,
      inputIdx=input_idx,
      outputIdx=0,
      control_indices=control_indices,
      control_volts=control_volts,
      minVoltage=-1.5,
      maxVoltage=1.5, 
      eq_steps=100_000,
      sim_steps=1_000_000,
      num_intervals=100,
      seed=np.random.randint(low=1, high=2**30),
      cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/config.txt"),
      acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/uniform_acceptors_0.txt"),
      don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/uniform_donors_0.txt"),
      ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/electrodes.txt"),
      save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse"),
      file_name=Path(f"bo_mse_data_type={data_type}_input_idx={input_idx-1}_func_type={func_type}_lambda_1=0.0_lambda_2=1.0_c_range=1.5"),
      BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
      SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
      OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )
  
  if data_type == 'data_vMB_ring_1e7':
    slurm_single_IV(
      numOfPoints=50,
      inputIdx=input_idx,
      outputIdx=0,
      control_indices=control_indices,
      control_volts=control_volts,
      minVoltage=-1.5,
      maxVoltage=1.5,
      eq_steps=100_000,
      sim_steps=1_000_000,
      num_intervals=100,
      seed=np.random.randint(low=1, high=2**30),
      cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
      acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt"),
      don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
      ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
      save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse"),
      file_name=Path(f"bo_mse_data_type={data_type}_input_idx={input_idx-1}_func_type={func_type}_lambda_1=0.0_lambda_2=1.0_c_range=1.5"),
      BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
      SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
      OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )
  
  if data_type == 'grad_towards_1e7':
    slurm_single_IV(
      numOfPoints=50,
      inputIdx=input_idx,
      outputIdx=0,
      control_indices=control_indices,
      control_volts=control_volts,
      minVoltage=-1.5,
      maxVoltage=1.5,
      eq_steps=100_000,
      sim_steps=1_000_000,
      num_intervals=100,
      seed=np.random.randint(low=1, high=2**30),
      cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/config.txt"),
      acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/vMB_gradient_towards_output.txt"),
      don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt"),
      ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/electrodes.txt"),
      save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse"),
      file_name=Path(f"bo_mse_data_type={data_type}_input_idx={input_idx-1}_func_type={func_type}_lambda_1=0.0_lambda_2=1.0_c_range=1.5"),
      BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
      SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
      OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
      )

  if data_type == 'grad_away_1e7':
    slurm_single_IV(
      numOfPoints=50,
      inputIdx=input_idx,
      outputIdx=0,
      control_indices=control_indices,
      control_volts=control_volts,
      minVoltage=-1.5,
      maxVoltage=1.5,
      eq_steps=100_000,
      sim_steps=1_000_000,
      num_intervals=100,
      seed=np.random.randint(low=1, high=2**30),
      cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/config.txt"),
      acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/vMB_gradient_away_output.txt"),
      don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/uniform_donors_1.txt"),
      ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/electrodes.txt"),
      save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_BO_mse"),
      file_name=Path(f"bo_mse_data_type={data_type}_input_idx={input_idx-1}_func_type={func_type}_lambda_1=0.0_lambda_2=1.0_c_range=1.5"),
      BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
      SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
      OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )
