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

def slurm_batch_from_single_state(
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
        "batchFromSingleState",
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
        fem_res,
        dist_type,
        eps,
        input_idx,
        output_idx,
        eq_steps,
        sim_steps,
        num_of_tasks,
        LHCSeed,
        threadBaseSeed,
        save_folder,
        file_name,
        ROOT,
        WS_DIR,
        BINARY,
        SH_SCRIPT
):

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
        f"--femRes={fem_res}",
        f"--distType={dist_type}",
        f"--eps={eps}",
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

    control_volts = [-0.9331,  1.0766,  0.1923, -0.2020,  0.9212,  0.2802]
    
    #control_volts = [1.0, 0, 0, -1.0, 0, 0]
    control_indices = [2, 3, 4, 5, 6, 7]
    seed=np.random.randint(low=0, high=2**30 - 1)
    ns = [0.5, 1.0, 2.0, 5.0, 10.0]#[0.5, 1.0, 1.5, 2.0, 3.0]

    slurm_single_IV(
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
        seed=seed,
        cfg="data/sm_batches_1e6_vonMises_beta/config.txt",
        acc_cfg=f"data/sm_batches_1e6_vonMises_beta/vonMises_beta_SM.txt",
        don_cfg="data/sm_batches_1e6_vonMises_beta/donors.txt",
        ele_cfg="data/sm_batches_1e6_vonMises_beta/electrodes.txt",
        save_folder="robustness_experiment_sigmoid",
        file_name=f"opt_sim",
        ROOT=ROOT,
        WS_DIR=WS_DIR,
        BINARY=BINARY,
        SH_SCRIPT=SH_SCRIPT
    )

    for _n in ns:
        slurm_single_IV(
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
            seed=seed,
            cfg="data/sm_batches_1e6_vonMises_beta/config.txt",
            acc_cfg=f"configs/vonMises_beta_SM_pert_ns={_n}.txt",
            don_cfg="data/sm_batches_1e6_vonMises_beta/donors.txt",
            ele_cfg="data/sm_batches_1e6_vonMises_beta/electrodes.txt",
            save_folder="robustness_experiment_sigmoid",
            file_name=f"opt_sim_pert={_n}",
            ROOT=ROOT,
            WS_DIR=WS_DIR,
            BINARY=BINARY,
            SH_SCRIPT=SH_SCRIPT
        )

    thetas = [0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]
    """ for theta in thetas:
        slurm_single_IV(
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
            seed=3123123,
            cfg="configs/config.txt",
            acc_cfg=f"configs/jiggled_acc_uniform_{theta}.txt",
            don_cfg="configs/don_uniform.txt",
            ele_cfg="configs/electrodes.txt",
            save_folder="example_curves",
            file_name=f"jiggled_acc_uniform_{theta}",
            ROOT=ROOT,
            WS_DIR=WS_DIR,
            BINARY=BINARY,
            SH_SCRIPT=SH_SCRIPT
        ) """

    num_batches = 200
    batch_size = 1_000
    """ for i in range(0, 150):
        time.sleep(0.1)
        slurm_batch_from_single_state(
            batch_size=batch_size,
            min_V=-1.5, max_V=1.5,
            input_idx=1, output_idx=0,
            eq_steps=10_000, sim_steps=1_000_000, num_of_tasks=100,
            LHCSeed=np.random.randint(low=0, high=2**20 - 1), threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
            cfg="configs/config.txt", acc_cfg="configs/vonMises_beta_SM.txt", don_cfg="configs/donors.txt", ele_cfg="configs/electrodes.txt",
            save_folder="sm_batches_1e6_vonMises_beta", file_name=f"batch_1e6_{i}",
            ROOT=ROOT, WS_DIR = WS_DIR, BINARY = BINARY, SH_SCRIPT = str(ROOT / "scripts" / "slurm" / "batch_script.sh")
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
        file_name="single_device_0",
        ROOT=ROOT,
        WS_DIR=WS_DIR,
        BINARY=BINARY,
        SH_SCRIPT=str(ROOT / "scripts" / "slurm" / "helix_single.sh")
    ) """

    """ slurm_batch_from_single_state(
        batch_size=500,
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

""" slurm_batch_of_independant_states(
    batch_size=100,
    min_V=-1.5,
    max_V=1.5,
    n_acceptors=200,
    n_electrodes=8,
    n_donors=3,
    radius=150.0,
    nu0=1.0,
    a=20.0,
    T=77.0,
    energy_disorder=0.01,
    electrode_width=60.0,
    min_hop_distance=1.5,
    max_hop_distance=60.0,
    fem_res=10000,
    dist_type="uniform",
    eps=0.6,
    input_idx=1,
    output_idx=0,
    eq_steps=10_000,
    sim_steps=1_000_000,
    num_of_tasks=100,
    LHCSeed=np.random.randint(low=0, high=2**20 - 1),
    threadBaseSeed=np.random.randint(low=0, high=2**20 - 1),
    save_folder="independant_state_batches",
    file_name="test_batch_of_independant_states",
    ROOT=ROOT,
    WS_DIR=WS_DIR,
    BINARY=BINARY,
    SH_SCRIPT=str(ROOT / "scripts" / "slurm" / "batch_script.sh")
) """