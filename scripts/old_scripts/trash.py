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

#best_theta = theta_samples[max_index]  # shape (5,)
best_theta = theta_vector.detach().cpu().numpy()
print("Best index:", max_index)
print("Best theta:", best_theta)
print("Best score:", scores[max_index])

# --- 2. Build the 4×7 input tensor for this theta ---
input_tensor = torch.zeros(len(gate_inputs), num_voltages)
control_voltages = torch.from_numpy(best_theta).float().expand(len(gate_inputs), -1)
for j, c_idx in enumerate(control_indices):
    input_tensor[:, c_idx] = control_voltages[:, j]

for k in range(len(gate_inputs)):
    input_tensor[k, input_index_0] = gate_inputs[k][0]
    input_tensor[k, input_index_1] = gate_inputs[k][1]
input_tensor = input_tensor.to(device)

model.eval()
with torch.no_grad():
    y = model(input_tensor).detach().cpu().squeeze().numpy()
I_00 = y[0]
I_01 = y[1]
I_10 = y[2]
I_11 = y[3]

vals = [I_00, I_01, I_10, I_11] 
print("Currents:", vals)
print(vals)
max_val = max(vals)
min_val = min(vals)
max_margin = max_val - min_val
y_min = min_val - 0.10 * max_margin
fig, ax = plt.subplots(2, 1)

for a in ax:
    a.set_box_aspect(1)
    a.tick_params(labelbottom=False)

l_w = 2.
plot_kwargs = {
    'lw': l_w,
    'color': 'k'
}
x = np.linspace(0, 1, len(vals) + 1)best_theta = theta_vector.detach().cpu().numpy()
ax[0].vlines(x=x[0], ymin=y_min, ymax=vals[0], **plot_kwargs)
ax[0].vlines(x=x[-1], ymin=y_min, ymax=vals[-1], **plot_kwargs)
for i in range(1, len(vals) + 1):
    ax[0].hlines(y=vals[i-1], xmin=x[i-1], xmax=x[i], **plot_kwargs)
for j in range(1, len(vals)):
    ax[0].vlines(x=x[j], ymin=vals[j-1], ymax=vals[j], **plot_kwargs)
ax[0].grid(alpha=.7)
ax[0].set_facecolor('lavender')

in_0 = [gate_inputs[0][0], gate_inputs[1][0], gate_inputs[2][0], gate_inputs[3][0]]
in_1 = [gate_inputs[0][1], gate_inputs[1][1], gate_inputs[2][1], gate_inputs[3][1]]

colors = ['r', 'green']
lws = [0.5, 2]
z_orders = [4, 3]
lss = ['dashed']

for l, in_volts in enumerate([in_0, in_1]):
    max_v = max(in_volts)
    min_v = min(in_volts)
    v_margin = max_v - min_v
    signal_min = min_v - 0.10 * v_margin
    ax[1].vlines(x=x[0], ymin=signal_min, ymax=in_volts[0], lw=lws[l], alpha=1, zorder=z_orders[l], ls='dashed', color=colors[l], label=rf'$V_{l}$')
    ax[1].vlines(x=x[-1], ymin=signal_min, ymax=in_volts[-1], lw=lws[l], alpha=1, zorder=z_orders[l], ls='dashed', color=colors[l])
    for i in range(1, len(vals) + 1):
        ax[1].hlines(y=in_volts[i-1], xmin=x[i-1], xmax=x[i], lw=lws[l], alpha=1, zorder=z_orders[l], ls='dashed', color=colors[l])
    for j in range(1, len(vals)):
        ax[1].vlines(x=x[j], ymin=in_volts[j-1], ymax=in_volts[j], lw=lws[l], alpha=1, zorder=z_orders[l], ls='dashed', color=colors[l])
    ax[1].grid(alpha=.7)
    ax[1].set_facecolor('lavender')
    ax[1].legend()

    import numpy as np
import matplotlib.pyplot as plt
import torch

# --- 1. Get best theta from the search ---
best_theta = theta_samples[max_index]  # shape (5,)

print("Best index:", max_index)
print("Best theta:", best_theta)

# --- 2. Build the 4×7 input tensor for this theta ---
input_tensor = torch.zeros(len(gate_inputs), num_voltages)

control_voltages = torch.from_numpy(best_theta).float().expand(len(gate_inputs), -1)
for j, c_idx in enumerate(control_indices):
    input_tensor[:, c_idx] = control_voltages[:, j]

for k in range(len(gate_inputs)):
    input_tensor[k, input_index_0] = gate_inputs[k][0]
    input_tensor[k, input_index_1] = gate_inputs[k][1]

# --- 3. Run the model to get currents for the 4 gate cases ---
model.eval()
with torch.no_grad():
    input_tensor = input_tensor.to(device)
    y = model(input_tensor).squeeze(1)        # shape (4,)
    y_cpu = y.cpu().numpy()                   # to numpy for plotting

print("Currents for best theta (A):", y_cpu)

# --- 4. Compute a suggested threshold & masks (optional, for visualization) ---
truth_tensor = torch.tensor(truth_outputs, device=device)

low_mask  = (truth_tensor == 0)
high_mask = (truth_tensor == 1)

low_curr  = y[low_mask].cpu().numpy()
high_curr = y[high_mask].cpu().numpy()

max_low  = low_curr.max()
min_high = high_curr.min()
margin   = min_high - max_low

print("max_low (OFF):", max_low)
print("min_high (ON):", min_high)
print("margin (A):", margin)

# If margin > 0, a natural threshold is halfway between them:
if margin > 0:
    Ith = 0.5 * (max_low + min_high)
else:
    Ith = None

# --- 5. Plot the logic gate behavior for this theta ---
x = np.arange(len(gate_inputs))
labels = [f"({gi[0]}, {gi[1]})" for gi in gate_inputs]

colors = ["tab:blue" if t == 0.0 else "tab:orange" for t in truth_outputs]

plt.figure(figsize=(6, 4))
plt.bar(x, y_cpu, color=colors, alpha=0.8)
plt.xticks(x, labels)
plt.xlabel("Logic inputs (V0, V1)")
plt.ylabel("Current (A)")
plt.title("Logic gate behavior for best θ")

# Compute a tight y-range that doesn't start at zero
y_min = y_cpu.min()
y_max = y_cpu.max()
pad = 0.05 * (y_max - y_min + 1e-12)  # 5% padding
plt.ylim(y_min - pad, y_max + pad)

# Color legend: OFF / ON
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="tab:blue", label="Target 0"),
    Patch(facecolor="tab:orange", label="Target 1"),
]

# Draw threshold line if we have a positive margin
if Ith is not None:
    plt.axhline(Ith, linestyle="--", color="k",
                label=f"Suggested threshold\nI_th ≈ {Ith:.2e} A")

plt.legend(handles=legend_elements, loc="best")
plt.tight_layout()
plt.show()