import numpy as np
import os
import time
import subprocess

from pathlib import Path

def slurm_single_IV(
        numOfPoints,
        numOfSamples,
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
        f"--numOfSamples={numOfSamples}",
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
        numOfPoints,
        numOfSamples,
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
        file_name
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
        "batchRun",
        f"--batchSize={batch_size}",
        f"--numOfPoints={numOfPoints}",
        f"--numOfSamples={numOfSamples}",
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

if __name__ == "__main__":

    ROOT = Path(__file__).resolve().parents[2]
    WS_DIR = Path("/gpfs/bwfor/work/ws/hd_gy283-my_data")
    print(ROOT)
    BINARY = ROOT / "build" / "kmc_project"
    SH_SCRIPT = ROOT / "scripts" / "slurm" / "single_curve.sh"

    #control_volts = [-100.5, 100.7, 110.2, -100.9, -110.0, 110.0]
    control_volts = [0, 0, 0, 0, 0, 0]
    control_indices = [2, 3, 4, 5, 6, 7]

    slurm_single_IV(
        numOfPoints=100,
        numOfSamples=10,
        inputIdx=0,
        outputIdx=1,
        control_indices=control_indices,
        control_volts=control_volts,
        minVoltage=-1.5,
        maxVoltage=1.5,
        eq_steps=10_000,
        sim_steps=100_000,
        num_intervals=100,
        seed=1200031,
        cfg="configs/config.txt",
        acc_cfg="configs/acceptors_1.txt",
        don_cfg="configs/donors.txt",
        ele_cfg="configs/electrodes.txt",
        save_folder="test_folder",
        file_name="test_file_32",
        ROOT=ROOT,
        WS_DIR=WS_DIR,
        BINARY=BINARY,
        SH_SCRIPT=SH_SCRIPT
    )