import numpy as np
import os
import time
import subprocess

from pathlib import Path

def create_dopant_configuration(radius, n_a, n_d, name_a, name_d, mode):

    acc_pos, don_pos = np.zeros((n_a, 2)), np.zeros((n_d, 2))
    angles_acc, angles_don = np.zeros(n_a), np.zeros(n_d)
    r_acc, r_don = np.zeros(n_a), np.zeros(n_d)

    if mode == "uniform":

        angles_acc = np.random.uniform(0, 2*np.pi, size=n_a)
        angles_don = np.random.uniform(0, 2*np.pi, size=n_d)
        r_acc = radius * np.sqrt(np.random.uniform(0, 1, size=n_a))
        r_don = radius * np.sqrt(np.random.uniform(0, 1, size=n_d))

        acc_pos[:, 0] = r_acc * np.cos(angles_acc)
        acc_pos[:, 1] = r_acc * np.sin(angles_acc)
        don_pos[:, 0] = r_don * np.cos(angles_don)
        don_pos[:, 1] = r_don * np.sin(angles_don)

    elif mode == "normal_truncated":

        mean = np.array([0.0, 0.0])
        cov = np.eye(2) * (radius / 3.0)**2

        i = 0
        while i < n_a:
            sample = np.random.multivariate_normal(mean=mean, cov=cov)
            r = np.linalg.norm(sample)
            if r <= radius:
                acc_pos[i] = sample
                r_acc[i] = r
                angles_acc[i] = np.arctan2(sample[1], sample[0])
                i += 1

        j = 0
        while j < n_d:
            sample = np.random.multivariate_normal(mean=mean, cov=cov)
            r = np.linalg.norm(sample)
            if r <= radius:
                don_pos[j] = sample
                r_don[j] = r
                angles_don[j] = np.arctan2(sample[1], sample[0])
                j += 1

    np.savetxt(f"configs/{name_a}.txt", acc_pos, fmt="%.6f", delimiter="\t")
    np.savetxt(f"configs/{name_d}.txt", don_pos, fmt="%.6f", delimiter="\t")

    return angles_acc, angles_don, r_acc, r_don

def jiggle_configuration(angles, radii, dtheta_max, radial_sigma, name):

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

    angles_acc, angles_don, r_acc, r_don = create_dopant_configuration(150.0, 200, 3, "acc_uniform", "don_uniform", "uniform")
    thetas = [0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]
    for theta in thetas:
        jiggle_configuration(angles_acc, r_acc, dtheta_max=theta, radial_sigma=0.0, name=f"jiggled_acc_uniform_{theta}")