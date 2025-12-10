import numpy as np
import time
import argparse
from pathlib import Path
from jobs import slurm_single_IV

parser = argparse.ArgumentParser()

parser.add_argument("--n_points", type=int)
parser.add_argument("--eq_steps", type=int)
parser.add_argument("--sim_steps", type=int)
parser.add_argument("--num_intervals", type=int)
parser.add_argument("--in_idx", type=int)
parser.add_argument("--out_idx", type=int)
parser.add_argument("--min_v", type=float)
parser.add_argument("--max_v", type=float)
parser.add_argument("--cfg", type=str)
parser.add_argument("--acc_cfg", type=str)
parser.add_argument("--don_cfg", type=str)
parser.add_argument("--ele_cfg", type=str)
parser.add_argument("--save_folder", type=str)
parser.add_argument("--file_name", type=str)
parser.add_argument("--c_indices", nargs='+')
parser.add_argument("--c_volts", nargs='+')

args = parser.parse_args()

if __name__ == '__main__':
    print(args.c_indices)
    print(args.c_volts)
    slurm_single_IV(
        numOfPoints=args.n_points,
        inputIdx=args.in_idx,
        outputIdx=args.out_idx,
        control_indices=args.c_indices,
        control_volts=args.c_volts,
        minVoltage=args.min_v,
        maxVoltage=args.max_v,
        eq_steps=args.eq_steps,
        sim_steps=args.sim_steps,
        num_intervals=args.num_intervals,
        seed=np.random.randint(0, 2 ** 31 - 1),
        cfg=Path(args.cfg),
        acc_cfg=Path(args.acc_cfg),
        don_cfg=Path(args.don_cfg),
        ele_cfg=Path(args.ele_cfg),
        save_folder=Path(args.save_folder),
        file_name=Path(args.file_name),
        BINARY=Path('/home/hd/hd_hd/hd_gy283/kmc_project/build/kmc_project'),
        SH_SCRIPT=Path('/home/hd/hd_hd/hd_gy283/kmc_project/scripts/slurm/single_curve.sh'),
        OUT_DIR=Path('/home/hd/hd_hd/hd_gy283/kmc_project/slurm_out')
    )
