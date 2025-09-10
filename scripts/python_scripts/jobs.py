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
        ROOT,
        WS_DIR,
        BINARY,
        SH_SCRIPT
):
    CFG_DIR = ROOT / f"{cfg}"
    ACC_DIR = ROOT / f"{acc_cfg}"
    DON_DIR = ROOT / f"{don_cfg}"
    ELE_DIR = ROOT / f"{ele_cfg}"
    DATA_DIR = Path(WS_DIR) / f"{save_folder}"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    SLURM_OUT_DIR = ROOT / "slurm_out"
    SLURM_OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file_name = SLURM_OUT_DIR / f"{file_name}.out"
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
        f"--cfg={str(CFG_DIR)}",
        f"--accCfg={str(ACC_DIR)}",
        f"--donCfg={str(DON_DIR)}",
        f"--eleCfg={str(ELE_DIR)}",
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

if __name__ == "__main__":

    ROOT = Path(__file__).resolve().parents[2]
    WS_DIR = Path("/gpfs/bwfor/work/ws/hd_gy283-my_data")
    print(ROOT)
    BINARY = ROOT / "build" / "kmc_project"
    SH_SCRIPT = ROOT / "scripts" / "slurm" / "single_curve.sh"

    control_volts = [-0.9331,  1.0766,  0.1923, -0.2020,  0.9212,  0.2802]
    
    #control_volts = [1.0, 0, 0, -1.0, 0, 0]
    control_indices = [2, 3, 4, 5, 6, 7]
    seed=np.random.randint(low=0, high=2**30 - 1)

    num_batches = 100
    batch_size = 1_000

    for i in range(100, 200):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            output_idx=0,
            eq_steps=10_000, sim_steps=10_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**20 - 1), threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
            cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/vMB_configs/config.txt"), 
            acc_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/vMB_configs/acceptors.txt"), 
            don_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/vMB_configs/donors.txt"), 
            ele_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/vMB_configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/data_vMB"), 
            file_name=f"batch_{i}",
            BINARY = Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"), 
            SH_SCRIPT =Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/batch_script.sh"), 
            SLURM_OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        )

    #c_indices = [2, 3, 4, 5, 6, 7]
    #c_volts = [ 0.2066,  0.6799, -0.1204, -0.6649,  0.1976,  0.9915]

    seed = np.random.randint(low=0, high=12341312312)
    nrands = [0, 3, 5, 10, 20, 50, 100]

    """ for nrand in nrands:
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
            seed=np.random.randint(low=1, high=142341),
            cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/config.txt"),
            acc_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/configs/robustness_to_acc/robustness_acc_nrand={nrand}.txt"),
            don_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/donors.txt"),
            ele_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/configs/electrodes.txt"),
            save_folder=Path("/gpfs/bwfor/work/ws/hd_gy283-my_data/robustness_acceptor_changes"),
            file_name=Path(f"acc_changed_nrand={nrand}"),
            BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
            SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
            OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
        ) """

    c_indices = [2, 3, 4, 5, 6, 7]
    c_volts = [-1.4962374834526107, -1.475426536296481, -0.947778084932544, 0.5938290599036549, -0.9337312498957281, 0.5894644757290364]
    """ slurm_single_IV(
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
        seed=np.random.randint(low=1, high=142341),
        cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/data/sm_batches_1e6_vonMises_beta/config.txt"),
        acc_cfg=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/data/sm_batches_1e6_vonMises_beta/acceptors.txt"),
        don_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/data/sm_batches_1e6_vonMises_beta/donors.txt"),
        ele_cfg=Path("/home/hd/hd_hd/hd_gy283/kmc_project/data/sm_batches_1e6_vonMises_beta/electrodes.txt"),
        save_folder=Path("/home/hd/hd_hd/hd_gy283/kmc_project"),
        file_name=Path(f"simulated_curve_vonMises_beta_no_noise"),
        BINARY=Path("/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project"),
        SH_SCRIPT=Path("/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh"),
        OUT_DIR=Path("/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out")
    ) """

    #noise = np.linspace(0.5, 10, 100)
    #for n in noise:
