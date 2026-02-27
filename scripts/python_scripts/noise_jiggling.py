import numpy as np
import os
import time
import subprocess
import pandas as pd

from pathlib import Path
from jobs import slurm_single_IV

if __name__ == '__main__':
    data_type = 'grad_towards_1e7'
    func_type = 'Sine'
    input_idx = 4
    maxl2 = 2.0
    seed = 42
    control_indices = [1, 2, 3, 5, 6, 7]

    top_K_num = 25
    noise_strengths = [1, 2, 5, 10, 15, 20]
    if data_type == 'uni_1e7':
        acceptor_path = np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data_recover')