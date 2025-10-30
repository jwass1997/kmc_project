import numpy as np
import os
import time
import subprocess

import matplotlib.pyplot as plt
import matplotlib.patches as ps
from visualize_hops import visualize_current, current_snapshot
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pathlib import Path

def create_config(n_a, n_d, n_e,
                  radius,
                  nu_0, a, T,
                  en_dis,
                  ele_width, max_h, min_h,
                  no_dim,
                  Nr, Nt,
                  name, save_dir):
    
    Path(save_dir).mkdir(exist_ok=True, parents=True)
    path_to_save = Path(save_dir) / f"{name}.txt"

    with path_to_save.open('w', encoding='utf-8', newline='\n') as f:
        f.write(f'nAcceptors {n_a}\n')
        f.write(f'nDonors {n_d}\n')
        f.write(f'nElectrodes {n_e}\n')
        f.write(f'radius {radius}\n')
        f.write(f'nu0 {nu_0}\n')
        f.write(f'a {a}\n')
        f.write(f'T {T}\n')
        f.write(f'energyDisorder {en_dis}\n')
        f.write(f'electrodeWidth {ele_width}\n')
        f.write(f'minHopDistance {min_h}\n')
        f.write(f'maxHopDistance {max_h}\n')
        f.write(f'noDimension {no_dim}\n')
        f.write(f'Nr {Nr}\n')
        f.write(f'Nt {Nt}\n')       



def create_dopant_configuration(radius, n_a, n_d, name_a, name_d, mode, eps, save_dir):

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

    elif mode == "mixed":

        assert 0.0 <= eps <= 1.0, "eps must be between 0 and 1"

        mean = np.array([0.0, 0.0])
        cov = np.eye(2) * (radius / 3.5) ** 2
        
        i = 0
        while i < n_a:
            u = np.random.uniform(0.0, 1.0)
            if u < (1-eps):
                r_acc[i] = radius*np.sqrt(np.random.uniform(0.0, 1.0))
                angles_acc[i] = np.random.uniform(0.0, 2*np.pi)

                

                acc_pos[i, 0] = r_acc[i]*np.cos(angles_acc[i])
                acc_pos[i, 1] = r_acc[i]*np.sin(angles_acc[i])
                i += 1
            else:
                sample = np.random.multivariate_normal(mean=mean, cov=cov)
                r = np.linalg.norm(sample)

                if r <= radius:
                    acc_pos[i] = sample
                    r_acc[i] = r
                    angles_acc[i] = np.arctan2(sample[1], sample[0])
                    i += 1
        j = 0
        while j < n_d:
            u = np.random.uniform(0.0, 1.0)
            if u < (1-eps):
                r_don[j] = radius*np.sqrt(np.random.uniform(0.0, 1.0))
                angles_don[j] = np.random.uniform(0.0, 2*np.pi)

                don_pos[j, 0] = r_don[j]*np.cos(angles_don[j])
                don_pos[j, 1] = r_don[j]*np.sin(angles_don[j])
                j += 1
            else:
                sample = np.random.multivariate_normal(mean=mean, cov=cov)
                r = np.linalg.norm(sample)

                if r <= radius:
                    don_pos[j] = sample
                    r_don[j] = r
                    angles_don[j] = np.arctan2(sample[1], sample[0])
                    j += 1
                                
    file_A = Path(save_dir) / f"{name_a}.txt"
    file_D = Path(save_dir) / f"{name_d}.txt"
    np.savetxt(f"{str(file_A)}", acc_pos, fmt="%.6f", delimiter="\t")
    np.savetxt(f"{str(file_D)}", don_pos, fmt="%.6f", delimiter="\t")

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

def triangle_from_src_to_dst(src, dst, base_width):
    """
    Return the 3 vertices of an isosceles triangle whose base is centered at src
    (length = base_width) and whose tip is at dst.
    """
    v = dst - src
    L = np.linalg.norm(v)
    if L == 0:
        return None  # same point, nothing to draw
    u = v / L                        # direction from src to dst
    n = np.array([-u[1], u[0]])      # unit normal (perpendicular)
    half = 0.5 * base_width

    p_left  = src + half * n
    p_right = src - half * n
    tip     = dst
    return np.vstack([p_left, p_right, tip])

def log_scaling(v, eps):

    positives = v > 0
    v_positives = v[positives]
    v_min = np.min(v_positives)
    v_max = np.max(v_positives)

    if eps is None:
        eps = 1e-6 * v_max

    a = np.log(v_min + eps)
    b = np.log(v_max + eps)
    denom = b - a

    x = np.where(v + eps > 0.0, np.log(v + eps), -np.inf)

    alpha = np.zeros_like(v, dtype=float)
    alpha[positives] = (np.log(v[positives] + eps) - a) / denom

    alpha = np.clip(alpha, 0.0, 1.0)
    
    return alpha

def current_distribution(ax, add_cbar, device_data, alpha_cutoff):
    
    """ Load data"""
    events = device_data['event_matrix']
    simulation_time = device_data['sim_time']
    acc_xy = device_data['acc_xy']
    don_xy = device_data['don_xy']
    ele_xy = device_data['ele_xy']
    energies = device_data['energies']

    n_A = acc_xy.shape[0]
    n_E = ele_xy.shape[0]

    net_events = (events - events.T)
    dir_rows, dir_cols = np.where(net_events > 0)
    #lower_half_currents = np.abs(np.tril(net_events, k=-1))
    abs_current = np.abs(net_events) / simulation_time

    min_current = abs_current.min()
    max_current = abs_current.max()
    #print(min_current)
    #print(max_current)
    #print(events.sum())

    abs_current_normalized = abs_current / max_current
    #abs_current_normalized = log_scaling(abs_current, None)

    acc_face_clr = (0, 0.6, 0, 1)

    """ Energies to color """
    cmap = plt.get_cmap('magma')
    norm = plt.Normalize(vmin=energies.min(), vmax=energies.max())
    colors = cmap(norm(energies))

    ax.set_aspect('equal', adjustable='box')

    """ Device boundary """
    radius = float(np.max(np.linalg.norm(acc_xy, axis=1)))
    device_boundary = ps.Circle(
        xy=(0.0, 0.0),
        radius=radius+0.15,
        fill=None,
        edgecolor='gray',
        zorder=4,
        lw=4
    )
    ax.add_patch(device_boundary)

    """ Electrodes """
    arcs = []
    th = np.linspace(0.0, 315, n_E)
    th_width = 12.5
    for i in range(n_E):
        """arcs.append(
            ps.Arc(
                xy=(0.0, 0.0),
                width=2*radius + 0.40,
                height=2*radius + 0.40,
                linewidth=8,
                color=colors[i + n_A],
                theta1=th[i]-th_width,
                theta2=th[i]+th_width,
                zorder=5
            )
        )"""
        arcs.append(
            ps.Wedge(
                center=(0, 0),
                r=radius + .4,
                fc=colors[i+n_A],
                ec='k',
                theta1=th[i]-th_width,
                theta2=th[i]+th_width,
                width = 0.4,
                zorder=5,
            )
        )
        
    for i, arc in enumerate(arcs):
        ax.add_patch(arc)
        #ax.text(x=(radius + 0.4)*np.cos(th[i]), y=(radius + 0.4)*np.sin(th[i]), s=rf'$V_{i}$')
        
    """ Acceptor and donors """
    node_radius = 0.13
    base_width = 2.0 * node_radius
    for i, pt in enumerate(acc_xy):
        ax.add_patch(
            ps.Circle(
                xy=(pt[0], pt[1]),
                radius=node_radius,
                zorder=4,
                color=colors[i]
            )
        )

    """ Plot currents between acceptors """
    for i in range(n_A):
        for j in range(n_A):
            a_val = abs_current_normalized[i, j]
            if a_val > alpha_cutoff:
                """ax.plot(
                    [acc_xy[i, 0], acc_xy[j, 0]],
                    [acc_xy[i, 1], acc_xy[j, 1]],
                    color='black',
                    alpha=abs_current_normalized[i, j],
                    zorder=3,
                    lw=1.5
                )"""

                if net_events[i, j] > 0:
                    src = acc_xy[i]
                    dest = acc_xy[j]
                    tri = triangle_from_src_to_dst(src, dest, base_width)
                    if tri is None:
                        continue
                    ax.add_patch(
                        ps.Polygon(
                            tri,
                            closed=True,
                            facecolor='black',
                            edgecolor='none',
                            alpha=a_val,#alpha=np.sqrt(a_val),
                            zorder=3
                        )
                    )
                elif net_events[i, j] < 0:
                    src = acc_xy[j]
                    dest = acc_xy[i]
                    tri = triangle_from_src_to_dst(src, dest, node_radius)
                    if tri is None:
                        continue
                    ax.add_patch(
                        ps.Polygon(
                            tri,
                            closed=True,
                            facecolor='black',
                            edgecolor='none',
                            alpha=a_val,#alpha=np.sqrt(a_val),
                            zorder=3
                        )
                    )

    cbar = None

    if add_cbar is not None:
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.3)
        fig = ax.figure
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(f'$E / k_b T$')
    
    ax.set_xlim(-radius - 0.5, radius + 0.5)
    ax.set_ylim(-radius - 0.5, radius + 0.5)
    ax.axis('off')
    
    return ax, cbar