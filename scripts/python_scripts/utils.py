import numpy as np
import os
import time
import subprocess

from pathlib import Path

def create_dopant_configuration(radius, n_a, n_d, name_a, name_d):
    
    acc_pos, don_pos = np.zeros(shape=(n_a, 2)), np.zeros(shape=(n_d, 2))

    angles_acc = np.random.uniform(low=0, high=2*np.pi, size=n_a)
    angles_don = np.random.uniform(low=0, high=2*np.pi, size=n_d)

    r_acc = radius * np.sqrt(np.random.uniform(low=0, high=1, size=n_a))
    r_don = radius * np.sqrt(np.random.uniform(low=0, high=1, size=n_d))

    with open(f"configs/{name_a}.txt", "w") as f:
        for i in range(n_a):
            acc_pos[i][0] = r_acc[i]*np.cos(angles_acc[i])
            acc_pos[i][1] = r_acc[i]*np.sin(angles_acc[i])
            f.write(f"{acc_pos[i][0]}\t{acc_pos[i][1]}\n")

    with open(f"configs/{name_d}.txt", "w") as f:
        for i in range(n_d):
            don_pos[i][0] = r_don[i]*np.cos(angles_don[i])
            don_pos[i][1] = r_don[i]*np.sin(angles_don[i])
            f.write(f"{don_pos[i][0]}\t{don_pos[i][1]}\n")

    return angles_acc, angles_don, r_acc, r_don

def jiggle_configuration(angles, radii, radius, dtheta_max, radial_sigma, name):

    N = angles.size
    new_positions = np.zeros(shape=(N, 2))

    delta_theta = np.random.uniform(low=-dtheta_max, high=dtheta_max, size=N)
    #delta_r = np.random.normal(loc=0.0, scale=radial_sigma, size=N)

    new_thetas = angles + delta_theta
    #new_radii = radii + delta_r

    new_positions[:, 0] = radii*np.cos(new_thetas)
    new_positions[:, 1] = radii*np.sin(new_thetas)

    with open(f"configs/{name}.txt", "w") as f:
        for i in range(N):
            f.write(f"{new_positions[i][0]}\t{new_positions[i][1]}\n")

if __name__ == "__main__":

    angles_acc, angles_don, r_acc, r_don = create_dopant_configuration(150.0, 200, 3, "acc_0", "don_0")
    jiggle_configuration(angles_acc, r_acc, radius=150.0, dtheta_max=0.25, radial_sigma=0.0, name="jiggled_acc_0")