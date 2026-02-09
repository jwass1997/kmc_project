import numpy as np
from pathlib import Path
from jobs import slurm_single_IV
from utils import create_dopant_configuration, create_config

num_geom = 20
num_electrodes = 8
electrode_pairs = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 0]
]

if __name__ == "__main__":
    """create_config(n_a=200, n_d=3, n_e=8,
                  radius=150.0,
                  nu_0=1.0, a=20.0, T=77.0,
                  en_dis=0.01,
                  ele_width=60.0, max_h=60.0, min_h=1.5,
                  no_dim=1,
                  Nr=257, Nt=1440,
                  name='test_config', save_dir='/home/hd/hd_hd/hd_gy283/kmc_project/configs')"""
    save_dir = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper')
    save_dir.mkdir(parents=True, exist_ok=True)

    #for i in range(num_geom):
    #    create_dopant_configuration(radius=150.0, n_a=200, n_d=3, name_a=f'acc_{i}', name_d=f'don_{i}', mode='uniform', eps=None, save_dir='/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper')
    en_dis = 0.0
    for j in electrode_pairs:
        input_idx, output_idx = j[1], j[0]
        c_indices = [c_idx for c_idx in range(num_electrodes) if c_idx not in (input_idx, output_idx)]
        for k in range(num_geom):
            slurm_single_IV(
                numOfPoints=50,
                inputIdx=input_idx,
                outputIdx=output_idx,
                control_indices=c_indices,
                control_volts=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                minVoltage=-1.5,
                maxVoltage=1.5,
                eq_steps=10_000,
                sim_steps=1_000_000,
                num_intervals=100,
                seed=np.random.randint(low=1, high=2**30),
                cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/configs/config_paper_ed={en_dis}.txt"),
                acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/configs/acc_{k}.txt"),
                don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/configs/don_{k}.txt"),
                ele_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt"),
                save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/en_dis={en_dis}"),
                file_name=Path(f"sim_input_idx={input_idx}_output_idx={output_idx}_geom={k}_en_dis={en_dis}"),
                BINARY=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                SH_SCRIPT=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                OUT_DIR=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
            )


    temperatures = np.linspace(77.0, 300.0, 10)

    """for temp in temperatures:
        create_config(n_a=200, n_d=3, n_e=8,
                  radius=150.0,
                  nu_0=1.0, a=20.0, T=temp,
                  en_dis=0.01,
                  ele_width=60.0, max_h=60.0, min_h=1.5,
                  no_dim=1,
                  Nr=257, Nt=1440,
                  name=f'cfg_{temp}', save_dir='/home/hd/hd_hd/hd_gy283/kmc_project/configs')
        resistance_data = Path(f'/gpfs/bwfor/work/ws/hd_gy283-my_data/resistance_from_paper/en_dis_{en_dis}/temp_{temp}')
        resistance_data.mkdir(exist_ok=True, parents=True)
        for j in electrode_pairs:
            input_idx, output_idx = j[1], j[0]
            c_indices = [c_idx for c_idx in range(num_electrodes) if c_idx not in (input_idx, output_idx)]
            for k in range(num_geom):
                slurm_single_IV(
                    numOfPoints=50,
                    inputIdx=input_idx,
                    outputIdx=output_idx,
                    control_indices=c_indices,
                    control_volts=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    minVoltage=-1.5,
                    maxVoltage=1.5,
                    eq_steps=10_000,
                    sim_steps=1_000_000,
                    num_intervals=100,
                    seed=np.random.randint(low=1, high=2**30),
                    cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/cfg_{temp}.txt"),
                    acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/configs/acc_{k}.txt"),
                    don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/iv_exp_from_paper/configs/don_{k}.txt"),
                    ele_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt"),
                    save_folder=resistance_data),
                    file_name=Path(f"sim_input_idx={input_idx}_output_idx={output_idx}_geom={k}_en_dis={en_dis}"),
                    BINARY=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
                    SH_SCRIPT=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
                    OUT_DIR=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
                )"""

    
        