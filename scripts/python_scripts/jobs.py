import numpy as np
import os
import time
import subprocess

from pathlib import Path

def slurm_single_device(
        control_indices,
        control_volts,
        output_idx,
        eq_steps,
        sim_steps,
        num_intervals,
        seed,
        cfg,
        acc_cfg,
        don_cfg,
        ele_cfg,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        OUT_DIR
):
    DATA_DIR = save_folder
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = Path(OUT_DIR) / f"{file_name}.out"
    params = {
        idx: c_v
        for idx, c_v in zip(control_indices, control_volts)
    }
    
    control_volts_args = []
    for idx, c_v in params.items():
        control_volts_args.append(f"--c_v={idx}={c_v}")

    args = [
        "sbatch",
        f"--output={output_file_name}",
        f"{str(SH_SCRIPT)}",
        f"{str(BINARY)}",
        f"singleRun",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numIntervals={num_intervals}",
        f"--seed={seed}",
        f"--cfg={str(cfg)}",
        f"--accCfg={str(acc_cfg)}",
        f"--donCfg={str(don_cfg)}",
        f"--eleCfg={str(ele_cfg)}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]
   
    cmd = args + control_volts_args

    result = subprocess.run(cmd, capture_output=True, text=True)
    #print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_single_IV(
        numOfPoints,
        inputIdx,
        outputIdx,
        control_indices,
        control_volts,
        minVoltage,
        maxVoltage,
        eq_steps,
        sim_steps,
        num_intervals,
        seed,
        cfg,
        acc_cfg,
        don_cfg,
        ele_cfg,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        OUT_DIR
):
    DATA_DIR = Path(save_folder)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = Path(OUT_DIR) / f"{file_name}.out"
    params = {
        idx: c_v
        for idx, c_v in zip(control_indices, control_volts)
    }
    
    control_volts_args = []
    for idx, c_v in params.items():
        control_volts_args.append(f"--c_v={idx}={c_v}")

    args = [
        "sbatch",
        f"--output={output_file_name}",
        f"{SH_SCRIPT}",
        f"{BINARY}",
        f"singleCurve",
        f"--numOfPoints={numOfPoints}",
        f"--inputIdx={inputIdx}",
        f"--outputIdx={outputIdx}",
        f"--minVoltage={minVoltage}",
        f"--maxVoltage={maxVoltage}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numIntervals={num_intervals}",
        f"--seed={seed}",
        f"--cfg={cfg}",
        f"--accCfg={acc_cfg}",
        f"--donCfg={don_cfg}",
        f"--eleCfg={ele_cfg}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]
   
    cmd = args + control_volts_args

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_batch_from_single_state(
        batch_size,
        min_V,
        max_V,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_tasks,
        LHCSeed,
        threadBaseSeed,
        cfg,
        acc_cfg,
        don_cfg,
        ele_cfg,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        SLURM_OUT_DIR
):
    DATA_DIR = Path(save_folder)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(SLURM_OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = SLURM_OUT_DIR / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batchFromSingleState",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfTasks={num_of_tasks}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--cfg={cfg}",
        f"--accCfg={acc_cfg}",
        f"--donCfg={don_cfg}",
        f"--eleCfg={ele_cfg}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_batch_from_single_state_sem(
        batch_size,
        min_V,
        max_V,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_devices,
        LHCSeed,
        threadBaseSeed,
        cfg,
        acc_cfg,
        don_cfg,
        ele_cfg,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        SLURM_OUT_DIR
):
    DATA_DIR = Path(save_folder)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(SLURM_OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = SLURM_OUT_DIR / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batchFromSingleState_sem",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfDevices={num_of_devices}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--cfg={cfg}",
        f"--accCfg={acc_cfg}",
        f"--donCfg={don_cfg}",
        f"--eleCfg={ele_cfg}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_batch_of_independant_states(
        batch_size,
        min_V,
        max_V,
        n_acceptors,
        n_electrodes,
        n_donors,
        radius,
        nu0,
        a,
        T,
        energy_disorder,
        electrode_width,
        min_hop_distance,
        max_hop_distance,
        Nr,
        Nt,
        dist_type,
        n_comps,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_tasks,
        LHCSeed,
        threadBaseSeed,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        SLURM_OUT_DIR
):

    DATA_DIR = Path(save_folder)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(SLURM_OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = Path(SLURM_OUT_DIR) / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batchOfMultipleStates",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--nAcceptors={n_acceptors}",
        f"--nElectrodes={n_electrodes}",
        f"--nDonors={n_donors}",
        f"--radius={radius}",
        f"--nu0={nu0}",
        f"--a={a}",
        f"--T={T}",
        f"--energyDisorder={energy_disorder}",
        f"--electrodeWidth={electrode_width}",
        f"--minHopDistance={min_hop_distance}",
        f"--maxHopDistance={max_hop_distance}",
        f"--Nr={Nr}",
        f"--Nt={Nt}",
        f"--distType={dist_type}",
        f"--n_comps={n_comps}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfTasks={num_of_tasks}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_batch_temp_dependent(
        batch_size,
        min_V,
        max_V,
        Tmin,
        Tmax,
        n_acceptors,
        n_electrodes,
        n_donors,
        radius,
        nu0,
        a,
        energy_disorder,
        electrode_width,
        min_hop_distance,
        max_hop_distance,
        Nr,
        Nt,
        dist_type,
        n_comps,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_tasks,
        LHCSeed,
        threadBaseSeed,
        save_folder,
        file_name,
        BINARY,
        SH_SCRIPT,
        SLURM_OUT_DIR
):

    DATA_DIR = Path(save_folder)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    Path(SLURM_OUT_DIR).mkdir(parents=True, exist_ok=True)
    output_file_name = Path(SLURM_OUT_DIR) / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batchOfSamplesWTemp",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--Tmin={Tmin}",
        f"--Tmax={Tmax}",
        f"--nAcceptors={n_acceptors}",
        f"--nElectrodes={n_electrodes}",
        f"--nDonors={n_donors}",
        f"--radius={radius}",
        f"--nu0={nu0}",
        f"--a={a}",
        f"--energyDisorder={energy_disorder}",
        f"--electrodeWidth={electrode_width}",
        f"--minHopDistance={min_hop_distance}",
        f"--maxHopDistance={max_hop_distance}",
        f"--Nr={Nr}",
        f"--Nt={Nt}",
        f"--distType={dist_type}",
        f"--n_comps={n_comps}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfTasks={num_of_tasks}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

if __name__ == "__main__":

    ROOT = Path(__file__).resolve().parents[2]
    WS_DIR = Path("/gpfs/bwfor/work/ws/hd_gy283-my_data")
    print(ROOT)
    BINARY = ROOT / "build" / "kmc_project"
    SH_SCRIPT = ROOT / "scripts" / "slurm" / "single_curve.sh"

    """for i in range(1, 2):
        time.sleep(0.1)
        slurm_batch_from_single_state_sem(
            batch_size=100,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=100_000, sim_steps=1_000_000, num_of_devices=10,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e6_with_sem_test_batch/configs/config.txt"), 
            acc_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e6_with_sem_test_batch/configs/uniform_acceptors_0.txt"), 
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e6_with_sem_test_batch/configs/uniform_donors_0.txt"), 
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e6_with_sem_test_batch/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e6_with_sem_test_batch"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"),   
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )"""

    # for i in range(200, 400):
    #     time.sleep(0.1)
    #     slurm_batch_of_independant_states(
    #         batch_size=500,
    #         min_V=-1.5,
    #         max_V=1.5,
    #         n_acceptors=200,
    #         n_electrodes=8,
    #         n_donors=3,
    #         radius=150.0,
    #         nu0=1.0,
    #         a=20.0,
    #         T=77.0,
    #         energy_disorder=0.01,
    #         electrode_width=60.0,
    #         min_hop_distance=1.5,
    #         max_hop_distance=60.0,
    #         Nr=257,
    #         Nt=1440,
    #         dist_type="uniform",
    #         n_comps=8,
    #         output_idx=0,
    #         eq_steps=100_000,
    #         sim_steps=10_000_000,
    #         num_of_tasks=100,
    #         LHCSeed=np.random.randint(low=0, high=2**31),
    #         threadBaseSeed=np.random.randint(low=0, high=2**31),
    #         save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/uni_cfg_samples_1e7"),
    #         file_name=f"batch_{i}",
    #         BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
    #         SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"),
    #         SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    #     )

    num_batches = 200
    batch_size = 1000
    """ for i in range(num_batches, 400):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=1_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/config.txt"), 
            acc_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/vMB_gradient_away_output.txt"), 
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/uniform_donors_1.txt"), 
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """
    
    """ for i in range(200):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=1_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**31 - 1), threadBaseSeed=np.random.randint(low=0, high=2**31 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/config.txt"), 
            acc_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/vMB_gradient_towards_output.txt"), 
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05/configs/uniform_donors_1.txt"), 
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/ddata_vMB_grad_away_sig=0.05/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_grad_away_sig=0.05"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """
    
    """ for i in range(200):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**30 - 1), threadBaseSeed=np.random.randint(low=0, high=2**30 - 1),
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_n_a=100/configs_n_a=100/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_n_a=100/configs_n_a=100/acceptors.txt"),
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_n_a=100/configs_n_a=100/donors.txt"),
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_n_a=100/configs_n_a=100/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_n_a=100"),
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """
    
    c_indices = [2, 3, 4, 5, 6, 7]
    #control_volts = np.random.uniform(low=-1.5, high=1.5, size=len(c_indices)).tolist()#[0.2, 1.1, -0.5, -1, 0.5, 0.7]
    control_volts = [0.6219740359908369, 1.5, 1.329525906935132, 1.5, 1.239992586471815, 0.7976947763969898]
    N = 2000

    steps = [100_000, 1_000_000]

    """for s in steps:
        slurm_single_IV(
            numOfPoints=100,
            inputIdx=1,
            outputIdx=0,
            control_indices=c_indices,
            control_volts=[0.3, -1.2, 0.2, -0.2, 0.5, .5, .8],
            minVoltage=0.3,
            maxVoltage=0.3,
            eq_steps=100_000,
            sim_steps=s,
            num_intervals=100,
            seed=np.random.randint(low=1, high=2**30),
            cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/config.txt"),
            acc_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/acceptors.txt"),
            don_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/donors.txt"),
            ele_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt"),
            save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets"),
            file_name=Path(f"steps={s}"),
            BINARY=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
            OUT_DIR=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )"""
    c_indices = [2, 3, 4, 5, 6, 7]
    c_volts = [-1.1192982948659904, 0.11237314251492486, -0.5453210281976734, -0.3224823080683447, 1.5, 0.860979202982826]
    """ slurm_single_IV(
        numOfPoints=50,
        inputIdx=1,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=control_volts,
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=10_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB_own_pde_solver/vMB_configs/config_0.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB_own_pde_solver/vMB_configs/acceptors.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB_own_pde_solver/vMB_configs/donors.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB_own_pde_solver/vMB_configs/electrodes.txt"),
        save_folder=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/relu"),
        file_name=Path(f"sim_0"),
        BINARY=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    ) """
    v = 0.0
    """for i in range(3):
        slurm_single_device(
            control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            control_volts = np.random.uniform(-1.5, 1.5, size=8).tolist(),
            output_idx=0,
            eq_steps=100_000,
            sim_steps=1_000_000,
            num_intervals=100,
            seed=np.random.randint(low=1, high=2**30),
            cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/uniform_acceptors_0.txt"),
            don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/uniform_donors_0.txt"),
            ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_a=5_ed=0.1_1e7/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets"),
            file_name=Path(f"device_{i}_a=5_ed=0.1"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/helix_single.sh"),
            OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )"""
    #c_volts = [-0.1381, -0.9105,  1.2631, -0.4615, -1.0557, -1.2425]
    """slurm_single_IV(
        numOfPoints=50,
        inputIdx=1,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=[ 0.6588,  0.6920,  0.9834, -1.0970,  0.3839,  0.6892],
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"vMB_sm_evaluation_torch_seed=44"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )"""

    """slurm_single_IV(
        numOfPoints=50,
        inputIdx=1,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=[ 1.1468,  1.2450, -0.3514,  1.3779, -0.3287,  0.3027],
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/vMB_gradient_away_output.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/uniform_donors_1.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_away_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"grad_away_sm_evaluation_torch_seed=42"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )"""

    """slurm_single_IV(
        numOfPoints=50,
        inputIdx=1,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=[ 0.6588,  0.6920,  0.9834, -1.0970,  0.3839,  0.6892],
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/uniform_acceptors_0.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/uniform_donors_0.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/uni_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"uni_sm_evaluation_44"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )"""
    # for i in [1000]:
    #     slurm_single_IV(
    #         numOfPoints=i,
    #         inputIdx=1,
    #         outputIdx=0,
    #         control_indices=c_indices,
    #         control_volts=[ 0.6588,  0.6920,  0.9834, -1.0970,  0.3839,  0.6892],
    #         minVoltage=-1.5,
    #         maxVoltage=1.5,
    #         eq_steps=100_000,
    #         sim_steps=10_000_000,
    #         num_intervals=100,
    #         seed=np.random.randint(low=1, high=2**30),
    #         cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/config.txt"),
    #         acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/vMB_gradient_towards_output.txt"),
    #         don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/uniform_donors_1.txt"),
    #         ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/grad_towards_1e7/configs/electrodes.txt"),
    #         save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
    #         file_name=Path(f"iv_speed_kmc_test_{i}"),
    #         BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
    #         SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
    #         OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    #     )

    """slurm_single_IV(
        numOfPoints=50,
        inputIdx=1,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=[ 0.0795349, 0.37506043, 1.06774043, -1.43881433, -1.25948799, -1.22106064],
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"data_vMB_ring_1e7_parabola"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )"""


    """ANNULAR"""
    data_type = 'data_vMB_ring_1e7'
    opt_type = 'random'
    func_type = 'linear'
    slurm_single_IV(
        numOfPoints=50,
        inputIdx=4,
        outputIdx=0,
        control_indices=c_indices,
        control_volts=[0.8753402530922632, -1.0253061239863164, -1.422887499677454, 0.3632564603904136, -0.5936841118833686, -0.4706212141182309],
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**30),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/found_1d_functions"),
        file_name=Path(f"opt_type={opt_type}_func_type={func_type}_data_type={data_type}"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )

    """ for n in noise_levels:
        seed = np.random.randint(low=1, high=2**15)
        slurm_single_IV(
            numOfPoints=100,
            inputIdx=1,
            outputIdx=0,
            control_indices=c_indices,
            control_volts=c_volts,
            minVoltage=-1.5,
            maxVoltage=1.5,
            eq_steps=10_000,
            sim_steps=1_000_000,
            num_intervals=100,
            seed=seed,
            cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/config.txt"),
            acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/perturbation_experiment_configs/acc_{n}.txt"),
            don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/donors.txt"),
            ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/perturbation_experiment"),
            file_name=Path(f"vMB_n={n}"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
            OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """
    volts = [0.0, -1.0, 0.2, 1.1, -0.5, -1, 0.5, 0.7]
    """slurm_single_device(
        control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        control_volts = [-0.19674096, -0.85058564, 0.0, -1.16885727, 1.5, -0.20346563,  1.46270051],
        output_idx=0,
        eq_steps=100_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=np.random.randint(low=1, high=2**31 - 1),
        cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/config.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/vMB_ring.txt"),
        don_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/uniform_donors_1.txt"),
        ele_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets/data_vMB_ring_1e7/configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"testing_current"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/helix_single.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    )"""

    """for steps in [100_000, 1_000_000, 10_000_000]:
        slurm_single_device(
            control_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            control_volts = volts,#np.random.uniform(low=-1.5, high=1.5, size=8).tolist(),
            eq_steps=100_000,
            sim_steps=f'{steps}',
            seed=np.random.randint(low=1, high=2**30),
            cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/config.txt"),
            acc_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/radial_distribution.txt"),
            don_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/donors_0.txt"),
            ele_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/final_datasets"),
            file_name=Path(f"steps_{steps}"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/helix_single.sh"),
            OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )"""