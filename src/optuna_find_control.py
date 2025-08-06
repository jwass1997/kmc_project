import numpy as np
import os
import subprocess
import optuna
import matplotlib.pyplot as plt

from pathlib import Path

def target_function(x: np.array) -> np.array:

    target = 1 / (1 + np.exp(4*x))

    return target

def mse(x: np.array, y:np.array) -> float:

    #N = x.shape[0]

    x_bar = np.mean(x)
    y_bar = np.mean(y)
    std_x = np.std(x, ddof=1)
    std_y = np.std(y, ddof=1)

    x_c = (x - x_bar) / std_x
    y_c = (y - y_bar) / std_y

    mse = np.mean((y_c - x_c)**2)

    return mse