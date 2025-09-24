import numpy as np
import os
import time
import subprocess

from pathlib import Path

def slurm_single_device(
        control_indices,
        control_volts,
        eq_steps,
        sim_steps,
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
    DATA_DIR = Path(WS_DIR) / f"{save_folder}"
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
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
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
    print("Running command:", " ".join(cmd))
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
        input_idx,
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
        "batchOfIndependantStates",
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
        f"--inputIdx={input_idx}",
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

    num_batches = 200
    batch_size = 1_000
    b_list = [73, 74, 113, 114, 115, 116]
    """ for i in b_list:#range(0, num_batches):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**30 - 1), threadBaseSeed=np.random.randint(low=0, high=2**30 - 1),
            cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/config.txt"), 
            acc_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/acceptors.txt"), 
            don_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/donors.txt"), 
            ele_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """
    
    c_indices = [2, 3, 4, 5, 6, 7]
    c_sig_volts = [0.8712577630625511, -0.5323602276247098, 0.9542597048359602, 1.3273684600791937, -0.3707830078391009, -0.726548944581742]
    c_sine_volts = [0.48431840576432994, 1.4355513682559609, 0.26524962777854855, 1.2232097957855352, 0.31940916409508846, 0.08077324262584651]
    control_volts = [-0.1, 1.1, -1.1, 0.1, -0.2, 1.1]
    slurm_single_IV(
        numOfPoints=100,
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
        cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/config_0.txt"),
        acc_cfg=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/acceptors.txt"),
        don_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/donors.txt"),
        ele_cfg=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB/vMB_configs/electrodes.txt"),
        save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/"),
        file_name=Path(f"vMB_sim_444"),
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