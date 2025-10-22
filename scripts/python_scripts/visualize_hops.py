import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches

from pathlib import Path
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

PYTHON_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = PYTHON_SCRIPT_DIR.parents[1]

def current_snapshot(hopping_counts, acceptor_pos, donor_pos, total_time):

    """
    Visualizes the net number of hops (current) between acceptor sites per unit time
    for a circular device. Line opacity (alpha) is mapped from current on a log scale.

    Parameters:
        hopping_counts : np.ndarray (n_acceptors x n_acceptors)
            Event count (hopping) matrix.
        acceptor_pos : np.ndarray of shape (n_acceptors, 2)
            (x,y) coordinates of acceptor sites.
        donor_pos : np.ndarray of shape (n_donors, 2)
            (x,y) coordinates of donor sites.
        total_time : float
            Total simulation time over which events were accumulated.
        alpha_min, alpha_max : float
            Min/max opacity for drawn edges.
        draw_colorbar : bool
            If True, show a grayscale colorbar with log ticks for |current|.
    """
    # --- Compute antisymmetric net hops per unit time
    net_hops = (hopping_counts - hopping_counts.T) / float(total_time)

    # --- Gather strictly positive magnitudes from lower triangle (exclude diagonal)
    lower = np.abs(np.tril(net_hops, k=-1)).ravel()
    pos = lower[lower > 0]

    # Early compute min/max among positive values (if any)
    if pos.size > 0:
        min_current = pos.min()
        max_current = pos.max()
    else:
        min_current = 0.0
        max_current = 0.0
    print(min_current, max_current)
    # --- Figure setup and device boundary
    max_radius = float(np.max(np.linalg.norm(acceptor_pos, axis=1)))
    padding = 0.1 * max_radius
    xlim = (-max_radius - padding, max_radius + padding)
    ylim = (-max_radius - padding, max_radius + padding)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    device_boundary = patches.Circle((0, 0), radius=max_radius + 0.15,
                                     fill=False, edgecolor="gray", zorder=4, lw=4)
    ax.add_patch(device_boundary)

    # --- Scatter donors and acceptors
    acc_face_color = (0, 0.6, 0, 0.8)
    acc_edge_color = (0, 0.6, 0, 1.0)
    if donor_pos is not None and len(donor_pos) > 0:
        ax.scatter(donor_pos[:, 0], donor_pos[:, 1], s=50, c="k", zorder=3, label="Donors")
    ax.scatter(acceptor_pos[:, 0], acceptor_pos[:, 1], s=50,
               facecolor=acc_face_color, edgecolors=acc_edge_color,
               zorder=3, label="Acceptors")

    alpha_min = 0.0
    alpha_max = 1.0
    # --- Draw edges with alpha mapped via log(current)
    n_acceptors = acceptor_pos.shape[0]
    for i in range(n_acceptors):
        for j in range(i + 1, n_acceptors):
            current_value = abs(net_hops[i, j])
            if current_value <= 0:
                continue  # zero (or negative before abs) -> no visible edge

            # Safe normalization in log domain:
            # handle degenerate cases (all positives equal, etc.)
            if min_current <= 0 or max_current <= 0:
                t = 1.0
            elif np.isclose(max_current, min_current):
                t = 1.0
            else:
                num = np.log(current_value) - np.log(min_current)
                den = np.log(max_current) - np.log(min_current)
                t = float(np.clip(num / den, 0.0, 1.0))

            alpha_val = float(alpha_min + t * (alpha_max - alpha_min))
            alpha_val = float(np.clip(alpha_val, alpha_min, alpha_max))  # final clamp

            ax.plot([acceptor_pos[i, 0], acceptor_pos[j, 0]],
                    [acceptor_pos[i, 1], acceptor_pos[j, 1]],
                    color="black", lw=1.5, alpha=alpha_val)

    # --- Legend
    ax.legend(loc="upper right")

    # --- Colorbar showing log10(|I|) on a linear scale (not log-scaled)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    if pos.size > 0 and (max_current > min_current) and (min_current > 0):
        vmin_log10 = np.log10(min_current)
        vmax_log10 = np.log10(max_current)
        norm_log10 = mcolors.Normalize(vmin=vmin_log10, vmax=vmax_log10)
        sm = ScalarMappable(norm=norm_log10, cmap=plt.cm.Greys)  # keeps your greyscale look
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(r"$\log_{10} |I_{ij}|$ (per unit time)")

        # integer tick marks if possible
        lo = np.ceil(vmin_log10)
        hi = np.floor(vmax_log10)
        if hi >= lo:
            ticks = np.arange(lo, hi + 1)
            cbar.set_ticks(ticks)

    elif pos.size > 0 and np.isclose(max_current, min_current) and (min_current > 0):
        # Degenerate case: all positives equal -> single log10 value
        val_log10 = np.log10(min_current)
        norm_log10 = mcolors.Normalize(vmin=val_log10 - 1e-9, vmax=val_log10 + 1e-9)
        sm = ScalarMappable(norm=norm_log10, cmap=plt.cm.Greys)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(rf"$\log_{{10}} |I_{{ij}}| \approx {val_log10:.3g}$")
        cbar.set_ticks([val_log10])

    else:
        # No positive currents; hide the colorbar axis
        cax.axis("off")

    plt.show()

def visualize_current(hopping_counts, acceptor_pos, donor_pos, total_time):
    """
    Visualizes the net number of hops (current) between sites per unit time
    for a circular device.
    
    Parameters:
        hopping_counts : np.ndarray
            The event count (hopping) matrix.
        acceptor_pos : np.ndarray of shape (n_acceptors, 2)
            The (x,y) coordinates of acceptor sites.
        donor_pos : np.ndarray of shape (n_donors, 2)
            The (x,y) coordinates of donor sites.
        electrode_pos : np.ndarray of shape (n_electrodes, 2)
            The (x,y) coordinates of electrodes.
        total_time : float
            Total simulation time over which events were accumulated.
    """
    # Compute net hops current (per unit time)
    net_hops = (hopping_counts - hopping_counts.T) / total_time

    # Get current magnitudes from lower triangular part of net_hops
    lower_half_vals = np.abs(np.tril(net_hops).flatten())
    min_current = lower_half_vals.min()
    max_current = lower_half_vals.max()

    # Create custom colormap for indicating intensity of net current
    current_range = np.linspace(0.0, max_current, 256)
    colors = np.array([[0, 0, 0, np.sqrt(val / max_current)] for val in current_range])
    custom_cmap = ListedColormap(colors)
    norm = mcolors.Normalize(vmin=min_current, vmax=max_current)
    sm = ScalarMappable(norm=norm, cmap=custom_cmap)
    sm.set_array([])

    # Define plot boundaries based on the acceptor positions (assume centered at (0, 0))
    max_radius = np.max(np.linalg.norm(acceptor_pos, axis=1))
    padding = 0.1 * max_radius
    xlim = (-max_radius - padding, max_radius + padding)
    ylim = (-max_radius - padding, max_radius + padding)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # Draw circular boundary representing the device outline
    device_boundary = patches.Circle((0, 0), radius=max_radius+0.15, fill=False, edgecolor="gray", zorder=4, lw=4)
    ax.add_patch(device_boundary)

    # Scatter plot of donor, acceptor, and electrode positions
    acc_face_color = (0, 0.6, 0, 0.8)
    acc_edge_color = (0, 0.6, 0, 1.0)
    ax.scatter(donor_pos[:, 0], donor_pos[:, 1], s=50, c="k", zorder=3, label="Donors")
    ax.scatter(acceptor_pos[:, 0], acceptor_pos[:, 1], s=50, facecolor=acc_face_color, edgecolors=acc_edge_color, zorder=3, label="Acceptors")

    # Create a colorbar to indicate net current intensity
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(sm, cax=cax)
    #cbar.set_label("")

    alpha_min = 0.1
    alpha_max = 1.0
    # Draw lines between each pair of acceptors with opacity based on current magnitude
    n_acceptors = acceptor_pos.shape[0]
    for i in range(n_acceptors):
        for j in range(i+1, n_acceptors):
            current_value = np.abs(net_hops[i, j])
            if current_value > 0:
                t = (np.log(current_value) - np.log(min_current)) / (np.log(max_current) - np.log(min_current) + 1e-12)
                t = np.clip(t, 0, 1)
                alpha_val = alpha_min + t * (alpha_max - alpha_min)
                ax.plot([acceptor_pos[i, 0], acceptor_pos[j, 0]],
                        [acceptor_pos[i, 1], acceptor_pos[j, 1]],
                        color="black", lw=1.5,
                        alpha=alpha_val)

    ax.legend(loc="upper right")
    #ax.set_title("Net Hops Current in Circular Device")
    plt.show()

if __name__ == "__main__":
    # Change the device name as needed.
    device = "6969"
    filename = f"/gpfs/bwfor/work/ws/hd_gy283-my_data/devices/single_device_0.npz"
    data = np.load(filename)

    # Load saved arrays
    hopping_counts = data["event_counts"]
    acceptor_coords = data["acceptor_coordinates"]
    donor_coords = data["donor_coordinates"]
    #electrode_coords = data["electrode_coordinates"] if "electrode_coordinates" in data else None
    total_time = data["device_time"]
    # total_time may be saved as an array or scalar; ensure we extract a scalar:
    if isinstance(total_time, np.ndarray):
        total_time = total_time.item()

    # Reshape coordinate arrays: they are stored flattened, so reshape to (-1,2)
    acceptor_pos = acceptor_coords.reshape(-1, 2)
    donor_pos = donor_coords.reshape(-1, 2)
    #electrode_pos = electrode_coords.reshape(-1, 2) if electrode_coords is not None else None

    # Visualize the net hops current on the circular device
    visualize_current(hopping_counts, acceptor_pos, donor_pos, total_time)

    # Optionally, to save the figure uncomment the line below and place it before plt.show() in visualizeCurrent
    plt.savefig(ROOT / "data" / "pngs" /f"single_device_0.png", dpi=300)
    #plt.show()