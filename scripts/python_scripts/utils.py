import torch
import numpy as np
import os
import time
import subprocess

import matplotlib.pyplot as plt
import matplotlib.patches as ps
from matplotlib.colors import LogNorm
from visualize_hops import visualize_current, current_snapshot
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import LogFormatterSciNotation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pathlib import Path

def normalize_over_range(x):
    
    z = (x - x.mean()) / (x.std() + 1e-8) 
    
    return z

def mse_loss_normalized(x, y):

    if isinstance(x, torch.Tensor):
        x = x.float()
    if isinstance(y, torch.Tensor):
        y = y.float()
        
    x_norm = (x - x.mean()) / (x.std() + 1e-8)
    y_norm = (y - y.mean()) / (y.std() + 1e-8)
    
    mse = torch.mean((y_norm - x_norm)**2).item()

    return mse

def f_sm(X_np: np.ndarray, model):
    with torch.no_grad():
        X = torch.tensor(X_np, dtype=torch.float32)
        y_full = model(X).squeeze(-1)
    return y_full[:, 0].detach().cpu().numpy()


def sample_X(N, min_v, max_v, d):
    return np.random.uniform(
        low=[min_v for i in range(d)],
        high=[max_v for i in range(d)],
        size=(N, d)
    )

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

def Beta_von_Mises_sample(alpha, beta, mu, kappa, radius, N):
    samples = []
    
    while len(samples) < N:
        # propose a batch of points (oversample a bit for efficiency)
        batch_size = max(100, N - len(samples))

        u = np.random.beta(a=alpha, b=beta, size=batch_size)
        r = radius * np.sqrt(u)
        theta = np.random.vonmises(mu=mu, kappa=kappa, size=batch_size)

        x = r * np.cos(theta)
        y = r * np.sin(theta)
        candidates = np.column_stack([x, y])
        
        # rejection step: keep only those inside the radius disc
        # (norm <= radius). This is equivalent to r <= radius.
        inside = np.linalg.norm(candidates, axis=1) <= radius
        accepted = candidates[inside]

        if accepted.size > 0:
            samples.append(accepted)

    # concatenate and trim to exactly N samples
    samples = np.vstack(samples)[:N]
    return samples
    
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

    e_0 = 1.60217663e-19
    current_scaling_factor = 1e9*1e12 # nanon ampere and time to secondes 

    net_events = (events - events.T)
    dir_rows, dir_cols = np.where(net_events > 0)
    #lower_half_currents = np.abs(np.tril(net_events, k=-1))
    abs_current = np.abs(net_events) / simulation_time * e_0 * current_scaling_factor

    positive = abs_current[abs_current > 0]
    vmin = positive.min()
    vmax = positive.max()
    curr_norm = LogNorm(vmin=vmin, vmax=vmax)
    curr_cmap = plt.get_cmap('Greys')
    
    def alpha_from_current(I, vmin, vmax, a_min=0.05, a_max=0.95):
        # map current to [0,1] in log space, then to [a_min, a_max]
        t = (np.log10(I) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
        t = np.clip(t, 0.0, 1.0)
        return a_min + t * (a_max - a_min)

    min_current = abs_current.min()
    max_current = abs_current.max()
    #print(min_current)
    #print(max_current)
    #print(events.sum())

    abs_current_normalized = abs_current #/ max_current
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
            I = abs_current[i, j]  # current in nA
            if I <= 0:
                continue
            
            # IMPORTANT: alpha_cutoff is now interpreted in nA
            if I < alpha_cutoff:
                continue
            
            a = alpha_from_current(I, vmin, vmax)
            fc = curr_cmap(curr_norm(I))  # color encodes current too (matches colorbar)
            
            if net_events[i, j] > 0:
                src = acc_xy[i]
                dest = acc_xy[j]
                tri = triangle_from_src_to_dst(src, dest, base_width)
            elif net_events[i, j] < 0:
                src = acc_xy[j]
                dest = acc_xy[i]
                tri = triangle_from_src_to_dst(src, dest, base_width)
            else:
                continue
            
            if tri is None:
                continue
            
            ax.add_patch(
                ps.Polygon(
                    tri,
                    closed=True,
                    facecolor=fc,      # was 'black'
                    edgecolor='none',
                    alpha=a,           # was a_val (way too large / not in [0,1])
                    zorder=3
                )
            )
            """a_val = abs_current_normalized[i, j]
            if a_val > alpha_cutoff:
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
                    )"""
    cbar_E = None
    cbar_I = None
    
    if add_cbar is not None:
        divider = make_axes_locatable(ax)
    
        # Put current colorbar on the RIGHT
        cax_I = divider.append_axes('right', size='5%', pad=0.25)
        sm_I = ScalarMappable(norm=curr_norm, cmap=curr_cmap)
        sm_I.set_array([])
        fig = ax.figure
        cbar_I = fig.colorbar(sm_I, cax=cax_I, format=LogFormatterSciNotation())
        cbar_I.set_label(r"Current (nA)")
    
        # Put energy colorbar further RIGHT (second bar)
        cax_E = divider.append_axes('right', size='5%', pad=0.90)
        sm_E = ScalarMappable(norm=norm, cmap=cmap)
        sm_E.set_array([])
        cbar_E = fig.colorbar(sm_E, cax=cax_E)
        cbar_E.set_label(r"$E / k_b T$")
        
    ax.set_xlim(-radius - 0.5, radius + 0.5)
    ax.set_ylim(-radius - 0.5, radius + 0.5)
    ax.axis('off')
    
    return ax, cbar_E, cbar_I

def pred_vs_true(ax, test_input, test_target, model, bins=140, use_log=True):

    test_input = test_input.to('cpu')
    y_true = test_target

    model = model.to('cpu')
    model.eval()
    with torch.no_grad():
        y_pred = model(test_input).detach().cpu()

    #y_pred = y_pred * model.y_std.numpy() + model.y_mean.numpy()

    y_pred = np.asarray(y_pred).ravel()
    y_true = np.asarray(y_true).ravel()

    pad = 0.1
    lo = min(y_pred.min(), y_true.min())
    hi = max(y_pred.max(), y_true.max())
    plot_range=[[lo, hi], [lo, hi]]

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    x_ticks = ax.get_xticks()
    ax.set_yticks(x_ticks)
    
    norm = LogNorm()

    #im = ax.hexbin(y_pred, y_true, bins=bins, gridsize=50, norm=norm, cmap="RdBu", zorder=5)
    _, _, _, im = ax.hist2d(y_pred, y_true, bins=bins, range=plot_range, norm=norm, cmap="RdBu", zorder=5)
    ax.plot([lo, hi], [lo, hi], linestyle='--', linewidth=2, color='k', alpha=0.7, zorder=6)

    fig = ax.figure
    cax = fig.add_axes([ax.get_position().x1 - 0.05, ax.get_position().y0, 0.04, ax.get_position().height])
    plt.colorbar(im, cax=cax)
    #ax.scatter(y_pred, y_true, s=4, edgecolors="blue", facecolor=None)
    ax.spines[:].set_linewidth(1.5)
    ax.tick_params(which='major', direction='in', width=1.5)
    ax.grid(linestyle='--', c='w')
    ax.set_facecolor('lightgray')

    ax.set_xlabel(r'$I_{\mathrm{pred}}$ $[ \nu_0 e_0 ]$')
    ax.set_ylabel(r'$I_{\mathrm{true}}$ $[ \nu_0 e_0 ]$')

    ax.set_aspect('equal', adjustable='box')

    return ax

def plot_gate_from_currents(gate_inputs, currents, currents_sem=None, sim_currents=None, sim_sem=None):
    y = np.asarray(currents, dtype=float).reshape(-1)
    if len(y) != len(gate_inputs):
        raise ValueError("`currents` length must match gate_inputs length.")

    if currents_sem is not None:
        sem = np.asarray(currents_sem, dtype=float).reshape(-1)
        if len(sem) != len(gate_inputs):
            raise ValueError("`currents_sem` length must match gate_inputs length.")
    else:
        sem = None

    if sim_currents is not None:
        y_sim = np.asarray(sim_currents, dtype=float).reshape(-1)
        if len(y_sim) != len(gate_inputs):
            raise ValueError("`sim_currents` length must match gate_inputs length.")
    else:
        y_sim = None

    if sim_sem is not None:
        sim_sem = 1.96*np.asarray(sim_sem, dtype=float).reshape(-1)
        if len(sim_sem) != len(gate_inputs):
            raise ValueError("`sim_sem` length must match gate_inputs length.")
    else:
        sim_sem = None

    in_0 = [gi[0] for gi in gate_inputs]
    in_1 = [gi[1] for gi in gate_inputs]

    fig, axes = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.12, top=0.98, hspace=0.15)

    labels = ["$V_0$", "$V_1$", "$I_{\\mathrm{out}} \\quad [nA]$"]
    t = np.arange(len(gate_inputs) + 1)  # 0..4

    def style_axis(ax, label):
        ax.set_ylabel(label, rotation=90.0)
        ax.grid(alpha=0.3)
        ax.set_facecolor("whitesmoke")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    # ---- output (mean step) ----
    y_out = np.r_[y, y[-1]]
    axes[2].step(t, y_out, where="post", lw=1, label="ground truth")

    # ---- shaded SEM band (mean ± sem) ----
    if sem is not None:
        lo = np.r_[y - sem, (y - sem)[-1]]
        hi = np.r_[y + sem, (y + sem)[-1]]
        axes[2].fill_between(t, lo, hi, step="post", alpha=0.25)

    # optional overlay: simulated
    if y_sim is not None:
        y_sim_out = np.r_[y_sim, y_sim[-1]]
        axes[2].step(t, y_sim_out, where="post", lw=2, ls="--", label="simulated")

        if sim_sem is not None:
            lo_s = np.r_[y_sim - sim_sem, (y_sim - sim_sem)[-1]]
            hi_s = np.r_[y_sim + sim_sem, (y_sim + sim_sem)[-1]]
            axes[2].fill_between(t, lo_s, hi_s, step="post", alpha=0.15)

        axes[2].legend(frameon=False, loc="best")

    style_axis(axes[2], labels[2])

    # ---- inputs ----
    y0 = np.r_[in_0, in_0[-1]]
    axes[0].step(t, y0, where="post", lw=2, c='r')
    style_axis(axes[0], labels[0])

    y1 = np.r_[in_1, in_1[-1]]
    axes[1].step(t, y1, where="post", lw=2, c='g')
    style_axis(axes[1], labels[1])

    centers = np.arange(len(gate_inputs)) + 0.5
    state_labels = ["00", "01", "10", "11"][:len(gate_inputs)]
    axes[-1].set_xticks(centers)
    axes[-1].set_xticklabels(state_labels)
    axes[-1].set_xlabel("Input $V_0, V_1 \\quad [V]$")

    return axes


"""def plot_gate_from_currents(gate_inputs, currents, sim_currents=None):

    y = np.asarray(currents, dtype=float).reshape(-1)
    if len(y) != len(gate_inputs):
        raise ValueError(f"`currents` length ({len(y)}) must match gate_inputs length ({len(gate_inputs)}).")

    if sim_currents is not None:
        y_sim = np.asarray(sim_currents, dtype=float).reshape(-1)
        if len(y_sim) != len(gate_inputs):
            raise ValueError(f"`sim_currents` length ({len(y_sim)}) must match gate_inputs length ({len(gate_inputs)}).")
    else:
        y_sim = None

    in_0 = [gi[0] for gi in gate_inputs]
    in_1 = [gi[1] for gi in gate_inputs]

    fig, axes = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.12, top=0.98, hspace=0.15)

    labels = ["$V_0$", "$V_1$", "$I_{\\mathrm{out}} \\quad [nA]$"]

    # time axis: 4 states -> 5 edges
    t = np.arange(len(gate_inputs) + 1)

    def style_axis(ax, label):
        ax.set_ylabel(label, rotation=90.0)
        ax.grid(alpha=0.3)
        ax.set_facecolor("whitesmoke")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    # output (ground truth)
    y_out = np.r_[y, y[-1]]
    axes[2].step(t, y_out, where="post", lw=2, label="ground truth")
    style_axis(axes[2], labels[2])

    # optional overlay: simulated currents
    if y_sim is not None:
        y_sim_out = np.r_[y_sim, y_sim[-1]]
        axes[2].step(t, y_sim_out, where="post", lw=2, ls="--", label="simulated")
        axes[2].legend(frameon=False, loc="best")

    # input 0
    y0 = np.r_[in_0, in_0[-1]]
    axes[0].step(t, y0, where="post", lw=2, c='r')
    style_axis(axes[0], labels[0])

    # input 1
    y1 = np.r_[in_1, in_1[-1]]
    axes[1].step(t, y1, where="post", lw=2, c='g')
    style_axis(axes[1], labels[1])

    # x-ticks at the middle of each interval, labeled by input combination
    centers = np.arange(len(gate_inputs)) + 0.5
    state_labels = ["00", "01", "10", "11"][:len(gate_inputs)]
    axes[-1].set_xticks(centers)
    axes[-1].set_xticklabels(state_labels)
    axes[-1].set_xlabel("Input $V_0, V_1 \\quad [V]$")

    # margins
    for ax in axes:
        ymin, ymax = ax.get_ylim()
        margin = 0.1 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - margin, ymax + margin)
        #ax.set_box_aspect(1)

    return axes"""

def plot_gate(model, theta, control_indices, gate_inputs):

    input_tensor = torch.zeros(len(gate_inputs), num_voltages)
    control_voltages = torch.from_numpy(theta).float().expand(len(gate_inputs), -1)
    for j, c_idx in enumerate(control_indices):
        input_tensor[:, c_idx] = control_voltages[:, j]
    
    for k in range(len(gate_inputs)):
        input_tensor[k, input_index_0] = gate_inputs[k][0]
        input_tensor[k, input_index_1] = gate_inputs[k][1]
    input_tensor = input_tensor.to(device)

    model.eval()
    with torch.no_grad():
        y = (model(input_tensor).detach().cpu().squeeze().numpy()).tolist()
    #print(y)
    in_0 = [gate_inputs[0][0], gate_inputs[1][0], gate_inputs[2][0], gate_inputs[3][0]]
    in_1 = [gate_inputs[0][1], gate_inputs[1][1], gate_inputs[2][1], gate_inputs[3][1]]
    
    fig, axes = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.12, top=0.98, hspace=0.15)

    labels = ["$V_0$", "$V_1$", "$I_{\mathrm{out}} \quad [nA]$"]
    
    # time axis: 4 states -> 5 edges
    t = np.arange(len(in_0) + 1)
    
    # helper: nice common style
    def style_axis(ax, label):
        ax.set_ylabel(label, rotation=90.0)#, ha="right", va="center")
        ax.grid(alpha=0.3)
        ax.set_facecolor("whitesmoke")
        # remove top/right spines
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    
    # output
    y_out = np.r_[y, y[-1]]
    axes[2].step(t, y_out, where="post", lw=2)
    style_axis(axes[2], labels[2])
    
    # input 0
    y0 = np.r_[in_0, in_0[-1]]
    axes[0].step(t, y0, where="post", lw=2, c='r')
    style_axis(axes[0], labels[0])
    
    # input 1
    y1 = np.r_[in_1, in_1[-1]]
    axes[1].step(t, y1, where="post", lw=2, c='g')
    style_axis(axes[1], labels[1])
    
    # x-ticks at the middle of each interval, labeled by input combination
    centers = np.arange(len(y)) + 0.5
    state_labels = ["00", "01", "10", "11"]
    axes[-1].set_xticks(centers)
    axes[-1].set_xticklabels(state_labels)
    axes[-1].set_xlabel("Input $V_0, V_1 \quad [V]$")
    
    # small margins so steps aren't glued to edges
    for ax in axes:
        ymin, ymax = ax.get_ylim()
        margin = 0.1 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - margin, ymax + margin)
        #ax.set_box_aspect(1)
    
    return axes