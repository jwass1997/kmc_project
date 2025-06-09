import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

def visualizeCurrent(hopping_counts, acceptor_pos, donor_pos, total_time):
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
    device_boundary = patches.Circle((0, 0), radius=max_radius, fill=False, edgecolor="gray", lw=2)
    ax.add_patch(device_boundary)

    # Scatter plot of donor, acceptor, and electrode positions
    ax.scatter(donor_pos[:, 0], donor_pos[:, 1], s=50, c="black", zorder=3, label="Donors")
    ax.scatter(acceptor_pos[:, 0], acceptor_pos[:, 1], s=50, c="green", zorder=3, label="Acceptors")

    # Create a colorbar to indicate net current intensity
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Net Current (|hops| per unit time)")

    # Draw lines between each pair of acceptors with opacity based on current magnitude
    n_acceptors = acceptor_pos.shape[0]
    for i in range(n_acceptors):
        for j in range(i+1, n_acceptors):
            current_value = np.abs(net_hops[i, j])
            if current_value > 0:
                ax.plot([acceptor_pos[i, 0], acceptor_pos[j, 0]],
                        [acceptor_pos[i, 1], acceptor_pos[j, 1]],
                        color="black", lw=1.5,
                        alpha=np.sqrt(current_value / max_current))

    ax.legend(loc="upper right")
    ax.set_title("Net Hops Current in Circular Device")
    plt.show()

if __name__ == "__main__":
    # Change the device name as needed.
    device = "device_201"
    filename = "../data/" + device + ".npz"
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
    visualizeCurrent(hopping_counts, acceptor_pos, donor_pos, total_time)

    # Optionally, to save the figure uncomment the line below and place it before plt.show() in visualizeCurrent
    plt.savefig(device + "_circular.png", dpi=300)
    #plt.show()