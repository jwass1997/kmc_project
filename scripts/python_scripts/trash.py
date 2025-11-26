from models import DeepSet

model_dict = torch.load(f'/gpfs/bwfor/work/ws/hd_gy283-my_data/models/sm_cond_mog_DS_nll.pth')

deep_set = DeepSet(in_dim=2, h_dim=90, num_outputs=1, out_dim=4, aggr_type='mean')
cond_sm = CondSM(layer_dims=[11, 90, 90, 90, 2], dropout_p=0.1, batch_norm=True)

deep_set.load_state_dict(model_dict['deep_set_state_dict'])
cond_sm.load_state_dict(model_dict['cond_sm_state_dict'])
print('model loaded')

X_pos = torch.from_numpy(np.loadtxt('/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni/uni_configs/acceptors.txt')).unsqueeze(0).float()
radius = 150.0
na = 200
R = np.sqrt(radius** 2 * np.pi / na )
X_pos = X_pos / radius

num_voltages = 7
num_outputs = 2

N = 100
voltage_range = torch.linspace(-1.5, 1.5, N)

v_quantiles = torch.arange(-1.5, 1.5, 0.3)

e_indices = [i for i in range(num_voltages)]

fig, ax = plt.subplots(nrows=2, ncols=num_voltages, figsize=(2.5 * num_voltages, 5), constrained_layout=True)

cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=v_quantiles[0], vmax=v_quantiles[-1])
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax.ravel().tolist(), label=rf'Voltage')

deep_set.eval()
cond_sm.eval()
for i in range(num_voltages):
    for v_q in v_quantiles:
        input_tensor = torch.empty(N, num_voltages, dtype=torch.float32)
        input_tensor[..., [col for col in range(num_voltages) if col != i]] = v_q
        input_tensor[:, i] = voltage_range
        z = deep_set(X_pos).squeeze(1).expand(N, -1)
        #T = torch.tensor(77.0, dtype=torch.float).reshape(1, 1).expand(N, -1)
        I = cond_sm(torch.cat([input_tensor, z], dim=1))
        mu = I[..., 0].detach().numpy()
        std = np.exp(I[..., 1].detach().numpy())

        if i == 0:
            ax[0, i].set_ylabel(r'$\mu$')
            ax[1, i].set_ylabel(r'$\sigma$')


        ax[0, i].plot(voltage_range, mu, color=cmap(norm(v_q)), lw=2.0)
        ax[1, i].plot(voltage_range, std, color=cmap(norm(v_q)), lw=2.0)

        ax[0, i].tick_params(labelbottom=False)
        ax[0, i].set_title(rf'$V_{i}$')
        
        ax[0, i].set_facecolor('lightgrey')
        ax[1, i].set_facecolor('lightgrey')

        ax[0, i].spines[:].set_linewidth(1.5)
        ax[1, i].spines[:].set_linewidth(1.5)

        ax[0, i].grid(True, which='both', c='w', ls='dashed')
        ax[1, i].grid(True, which='both', c='w', ls='dashed')

        ax[1, 0].set_xlabel(r'$V$') 