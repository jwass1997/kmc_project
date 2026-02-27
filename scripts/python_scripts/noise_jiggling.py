import numpy as np
import os
import time
import subprocess
import pandas as pd

from pathlib import Path
from jobs import slurm_single_IV

def add_noise_project_to_disk(xy: np.ndarray, sigma: float, R: float, seed: int | None = None) -> np.ndarray:
    """
    Add isotropic Gaussian noise N(0, sigma^2) to each (x,y), then project any points outside
    the disk of radius R back onto the boundary (same angle, radius clipped to R).

    xy: (N,2) array
    sigma: noise std dev in same units as xy
    R: disk radius
    """
    rng = np.random.default_rng(seed)
    noisy = xy + rng.normal(0.0, sigma, size=xy.shape)

    r = np.linalg.norm(noisy, axis=1)
    outside = r > R
    if np.any(outside):
        # scale radii down to exactly R, keep angles
        noisy[outside] *= (R / r[outside])[:, None]
    return noisy


def add_noise_resample_in_disk(xy: np.ndarray, sigma: float, R: float, seed: int | None = None,
                               max_tries: int = 10_000) -> np.ndarray:
    """
    Add isotropic Gaussian noise but *reject* any jiggled positions that fall outside the disk,
    resampling noise for those points until they are inside.

    This preserves the intended noise distribution better than projection near the boundary.
    """
    rng = np.random.default_rng(seed)
    out = xy.copy()

    remaining = np.ones(len(xy), dtype=bool)
    tries = 0

    while np.any(remaining):
        tries += 1
        if tries > max_tries:
            raise RuntimeError(f"Exceeded max_tries={max_tries}. "
                               "Try smaller sigma or use projection method.")
        # propose new positions for remaining points
        proposal = xy[remaining] + rng.normal(0.0, sigma, size=(remaining.sum(), 2))
        r = np.linalg.norm(proposal, axis=1)
        ok = r <= R

        # write accepted proposals into output
        idx = np.flatnonzero(remaining)
        out[idx[ok]] = proposal[ok]

        # keep trying for those still outside
        remaining[idx[ok]] = False

    return out

if __name__ == '__main__':

    data_type = 'data_vMB_ring_1e7'
    func_type = 'Parabola'
    maxl2 = 2.0
    seed = 42
    input_idx = 4
    
    control_indices = [1, 2, 3, 5, 6, 7]
    R = 150.0

    top_K_num = 5
    noise_strengths = [2, 5, 10, 20]
    num_noise_samples = 10
    base_seed = 42

    if data_type == 'uni_1e7':
        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )
        ref_dopants = np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/uniform_acceptors_0.txt')
        for sig in noise_strengths:
            save_path = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{data_type}/configs/jiggled_configs/sigma={sig}")
            save_path.mkdir(exist_ok=True, parents=True)
            for j in range(num_noise_samples):
                dopants_noisy = add_noise_project_to_disk(
                    ref_dopants, 
                    sigma=sig,
                    R=R, 
                    seed=base_seed + 1000*noise_strengths.index(sig) + j
                    )
                current_dopant_path = save_path / f"sig={sig}_jiggled_{j}.txt"
                np.savetxt(current_dopant_path, dopants_noisy)
                for k in range(top_K_num, 10):
                    job_seed = 10_000_000 + 1_000*noise_strengths.index(sig) + 10*j + k
                    time.sleep(0.2)
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
                        seed=job_seed,
                        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/config.txt"),
                        acc_cfg=Path(current_dopant_path),
                        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/uniform_donors_0.txt"),
                        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/uni_1e7/configs/electrodes.txt"),
                        save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}_sig={sig}"),
                        file_name=Path(f"topk={k}_jiggled={j}"),
                        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                    )
    if data_type == 'data_vMB_ring_1e7':
        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )
        ref_dopants = np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt')
        for sig in noise_strengths:
            save_path = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{data_type}/configs/jiggled_configs/sigma={sig}")
            save_path.mkdir(exist_ok=True, parents=True)
            for j in range(num_noise_samples):
                dopants_noisy = add_noise_project_to_disk(
                    ref_dopants, 
                    sigma=sig, 
                    R=R, 
                    seed=base_seed + 1000*noise_strengths.index(sig) + j
                    )
                current_dopant_path = save_path / f"sig={sig}_jiggled_{j}.txt"
                np.savetxt(current_dopant_path, dopants_noisy)
                for k in range(top_K_num, 10):
                    job_seed = 10_000_000 + 1_000*noise_strengths.index(sig) + 10*j + k
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
                        seed=job_seed,
                        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
                        acc_cfg=Path(current_dopant_path),
                        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
                        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
                        save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}_sig={sig}"),
                        file_name=Path(f"topk={k}_jiggled={j}"),
                        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                    )

    if data_type == 'grad_towards_1e7':
        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )
        ref_dopants = np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/vMB_gradient_towards_output.txt')
        for sig in noise_strengths:
            save_path = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{data_type}/configs/jiggled_configs/sigma={sig}")
            save_path.mkdir(exist_ok=True, parents=True)
            for j in range(num_noise_samples, 20):
                dopants_noisy = add_noise_project_to_disk(
                    ref_dopants, 
                    sigma=sig, 
                    R=R, 
                    seed=base_seed + 1000*noise_strengths.index(sig) + j
                    )
                current_dopant_path = save_path / f"sig={sig}_jiggled_{j}.txt"
                np.savetxt(current_dopant_path, dopants_noisy)
                for k in range(top_K_num):
                    job_seed = 10_000_000 + 1_000*noise_strengths.index(sig) + 10*j + k
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
                        seed=job_seed,
                        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/config.txt"),
                        acc_cfg=Path(current_dopant_path),
                        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt"),
                        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_towards_1e7/configs/electrodes.txt"),
                        save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}_sig={sig}"),
                        file_name=Path(f"topk={k}_jiggled={j}"),
                        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                    )

    if data_type == 'grad_away_1e7':
        df = np.loadtxt(
            f'/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/'
            f'topK_random_affine_data_type={data_type}_func={func_type}_cmin=-1.5_cmax=1.5_seed={seed}_maxl2={maxl2}_input_idx=3_K=100.txt'
            )
        ref_dopants = np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/vMB_gradient_away_output.txt')
        for sig in noise_strengths:
            save_path = Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/{data_type}/configs/jiggled_configs/sigma={sig}")
            save_path.mkdir(exist_ok=True, parents=True)
            for j in range(num_noise_samples, 20):
                dopants_noisy = add_noise_project_to_disk(
                    ref_dopants, 
                    sigma=sig, 
                    R=R, 
                    seed=base_seed + 1000*noise_strengths.index(sig) + j
                    )
                current_dopant_path = save_path / f"sig={sig}_jiggled_{j}.txt"
                np.savetxt(current_dopant_path, dopants_noisy)
                for k in range(top_K_num):
                    job_seed = 10_000_000 + 1_000*noise_strengths.index(sig) + 10*j + k
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
                        seed=job_seed,
                        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/config.txt"),
                        acc_cfg=Path(current_dopant_path),
                        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/uniform_donors_1.txt"),
                        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/final_datasets/grad_away_1e7/configs/electrodes.txt"),
                        save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data_recover/1d_funcs_random_mse/jiggled_topks_{data_type}_{func_type}_sig={sig}"),
                        file_name=Path(f"topk={k}_jiggled={j}"),
                        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                    )
