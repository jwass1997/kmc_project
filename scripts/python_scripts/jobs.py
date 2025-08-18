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

def slurm_single_batch(
        batch_size,
        min_V,
        max_V,
        input_idx,
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
        ROOT,
        WS_DIR,
        BINARY,
        SH_SCRIPT
):
    CFG_DIR = Path(ROOT) / f"{cfg}"
    ACC_DIR = Path(ROOT) / f"{acc_cfg}"
    DON_DIR = Path(ROOT) / f"{don_cfg}"
    ELE_DIR = Path(ROOT) / f"{ele_cfg}"
    DATA_DIR = Path(WS_DIR) / f"{save_folder}"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    SLURM_OUT_DIR = Path(ROOT) / "slurm_out"
    SLURM_OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file_name = SLURM_OUT_DIR / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batch",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--inputIdx={input_idx}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfTasks={num_of_tasks}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--cfg={str(CFG_DIR)}",
        f"--accCfg={str(ACC_DIR)}",
        f"--donCfg={str(DON_DIR)}",
        f"--eleCfg={str(ELE_DIR)}",
        f"--saveFolder={str(DATA_DIR)}",
        f"--fileName={file_name}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Running command:", " ".join(cmd))
    if result.returncode != 0:
        print("sbatch failed with stderr:\n", result.stderr)
    else:
        print("sbatch submission output:\n", result.stdout)

def slurm_single_batch_with_dist_param(
        batch_size,
        min_V,
        max_V,
        input_idx,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_tasks,
        LHCSeed,
        threadBaseSeed,
        dist_type,
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
    CFG_DIR = Path(ROOT) / f"{cfg}"
    ACC_DIR = Path(ROOT) / f"{acc_cfg}"
    DON_DIR = Path(ROOT) / f"{don_cfg}"
    ELE_DIR = Path(ROOT) / f"{ele_cfg}"
    DATA_DIR = Path(WS_DIR) / f"{save_folder}"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    SLURM_OUT_DIR = Path(ROOT) / "slurm_out"
    SLURM_OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file_name = SLURM_OUT_DIR / f"{file_name}.out"

    cmd = [
        "sbatch",
        f"--output={output_file_name}",
        str(SH_SCRIPT),
        str(BINARY),
        "batch_with_dist_param",
        f"--batchSize={batch_size}",
        f"--minVoltage={min_V}",
        f"--maxVoltage={max_V}",
        f"--inputIdx={input_idx}",
        f"--outputIdx={output_idx}",
        f"--eqSteps={eq_steps}",
        f"--simSteps={sim_steps}",
        f"--numOfTasks={num_of_tasks}",
        f"--LHCSeed={LHCSeed}",
        f"--threadBaseSeed={threadBaseSeed}",
        f"--distType={dist_type}",
        f"--cfg={str(CFG_DIR)}",
        f"--accCfg={str(ACC_DIR)}",
        f"--donCfg={str(DON_DIR)}",
        f"--eleCfg={str(ELE_DIR)}",
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

    control_volts = [1.4268, -0.5445, -1.1118, -1.4234, -0.4153,  1.4797]
    #control_volts = [1.0, 0, 0, -1.0, 0, 0]
    control_indices = [2, 3, 4, 5, 6, 7]

    """ slurm_single_IV(
        numOfPoints=100,
        inputIdx=1,
        outputIdx=0,
        control_indices=control_indices,
        control_volts=control_volts,
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=10_000,
        sim_steps=1_000_000,
        num_intervals=100,
        seed=4214124,
        cfg="configs/config.txt",
        acc_cfg="configs/vonMises_beta_SM.txt",
        don_cfg="configs/donors.txt",
        ele_cfg="configs/electrodes.txt",
        save_folder="example_curves",
        file_name="test_sine",
        ROOT=ROOT,
        WS_DIR=WS_DIR,
        BINARY=BINARY,
        SH_SCRIPT=SH_SCRIPT
    ) """

    num_batches = 200
    batch_size = 1_000
    """ for i in range(0, 150):
        time.sleep(0.1)
        slurm_single_batch(
            batch_size=batch_size,
            min_V=-1.5,
            max_V=1.5,
            input_idx=1,
            output_idx=0,
            eq_steps=10_000,
            sim_steps=1_000_000,
            num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**20 - 1),
            threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
            cfg="configs/config.txt",
            acc_cfg="configs/vonMises_beta_SM.txt",
            don_cfg="configs/donors.txt",
            ele_cfg="configs/electrodes.txt",
            save_folder="sm_batches_1e6_vonMises_beta",
            file_name=f"batch_1e6_{i}",
            ROOT=ROOT,
            WS_DIR = WS_DIR,
            BINARY = BINARY,
            SH_SCRIPT = str(ROOT / "scripts" / "slurm" / "batch_script.sh")
        ) """
    voltages = [0.1, -0.2, 1.2, 0.8, 0.6, -0.9, 0.5, 0.5]
    voltage_indices = [0, 1, 2, 3, 4, 5, 6, 7]
    """ slurm_single_device(
        control_indices=voltage_indices,
        control_volts=voltages,
        eq_steps=10_000,
        sim_steps=1_000_000,
        seed=np.random.randint(low=0, high=2**28),
        cfg="configs/config.txt",
        acc_cfg="configs/acceptors.txt",
        don_cfg="configs/donors.txt",
        ele_cfg="configs/electrodes.txt",
        save_folder="devices",
        file_name="single_device_eps=0.1",
        ROOT=ROOT,
        WS_DIR=WS_DIR,
        BINARY=BINARY,
        SH_SCRIPT=str(ROOT / "scripts" / "slurm" / "helix_single.sh")
    ) """

    slurm_single_batch_with_dist_param(
        batch_size=100,
        min_V=-1.5,
        max_V=1.5,
        input_idx=1,
        output_idx=0,
        eq_steps=10_000,
        sim_steps=1_000_000,
        num_of_tasks=100,
        LHCSeed=np.random.randint(low=0, high=2**20 - 1),
        threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
        dist_type="mixed",
        cfg="configs/config.txt",
        acc_cfg="configs/acceptors.txt",
        don_cfg="configs/donors.txt",
        ele_cfg="configs/electrodes.txt",
        save_folder="batches_with_dist_param",
        file_name=f"batch_1e6_dist_param",
        ROOT=ROOT,
        WS_DIR = WS_DIR,
        BINARY = BINARY,
        SH_SCRIPT = str(ROOT / "scripts" / "slurm" / "batch_script.sh")
    )

    """ slurm_single_batch(
        batch_size=100,
        min_V=-1.5,
        max_V=1.5,
        input_idx=1,
        output_idx=0,
        eq_steps=10_000,
        sim_steps=1_000_000,
        num_of_tasks=100,
        LHCSeed=np.random.randint(low=0, high=2**20 - 1),
        threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
        cfg="configs/config.txt",
        acc_cfg="configs/acceptors.txt",
        don_cfg="configs/donors.txt",
        ele_cfg="configs/electrodes.txt",
        save_folder="test_batches",
        file_name=f"test_batch",
        ROOT=ROOT,
        WS_DIR = WS_DIR,
        BINARY = BINARY,
        SH_SCRIPT = str(ROOT / "scripts" / "slurm" / "batch_script.sh")
    ) """