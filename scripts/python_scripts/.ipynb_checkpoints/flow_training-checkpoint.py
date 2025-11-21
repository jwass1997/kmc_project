import torch
import torch.nn as nn
import numpy as np
import normflows as nf

from tqdm import tqdm
from pathlib import Path
from torch.utils.data import Dataset, TensorDataset, DataLoader
from torch.utils.data import random_split
from train_sm import MakeDataset

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

num_batches = 200

data_dir = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_samples_1e6/')

coord_list = []
current_list = []
voltage_list = []
input_idx = 0 # input_idx AFTER sorting out the output electrode
output_idx = 0

if __name__ == "__main__":

    # --- Load data --- #
    
    for i in range(num_batches):
    
        fp = data_dir / f'batch_{i}.npz'
        try:
            batch = np.load(fp)
        except FileNotFoundError:
            continue        
    
        voltage_list.append(batch['inputs'])
        coord_list.append(batch['acc_xy'])    
    
        if batch['currents'].ndim != 1:
            current_list.append(batch['currents'][:, output_idx])
        else:
            current_list.append(batch['currents'])
    
    # You need to change this if you changed simulation parameters
    radius = 150.0
    n_a = 200
    R = np.sqrt(radius ** 2 * np.pi / n_a)
    print(R)
    
    coord_set = np.concatenate(coord_list) * R / radius
    current_set = np.concatenate(current_list)
    voltage_set = np.concatenate(voltage_list)
    
    print('Dataset loaded')
    
    print(len(coord_set))
    print(len(voltage_set))
    print(len(current_set))

    # --- Dataloaders --- #
        
    train_size = 0.8
    val_size = 0.1
    test_size = 0.1
    
    r2 = np.sum(coord_set**2, axis=-1)
    print(r2.max())
    idx = np.argsort(r2, axis=1)
    sorted_coord_set = np.take_along_axis(
        coord_set,
        idx[..., None],
        axis=1
    )
    
    X = torch.from_numpy(sorted_coord_set).float()
    C = torch.from_numpy(np.concatenate([voltage_set, current_set.reshape(-1, 1)], axis=1)).float()
    
    X = X[:200000, ...]
    C = C[:200000, ...]
    
    print(current_set.max(), current_set.min())
    
    print(X.shape)
    print(C.shape)
    
    dataset = TensorDataset(X, C)
    
    generator_1 = torch.Generator().manual_seed(32)
    train_ds, val_ds, test_ds = random_split(dataset, lengths=[train_size, val_size, test_size], generator=generator_1)
    
    batch_size = 2048
    num_workers = 2
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=num_workers, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)
    test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)
    
    print(len(train_loader.dataset))
    print(len(val_loader.dataset))
    print(len(test_loader.dataset))
    print('Dataloaders ready')

    # --- Flow model --- #
    
    K = 16

    latent_size = 400
    hidden_units = 128
    hidden_layers = 2
    
    flows = []
    for i in range(K):
        flows += [nf.flows.AutoregressiveRationalQuadraticSpline(latent_size, hidden_layers, hidden_units)]
        flows += [nf.flows.LULinearPermute(latent_size)]
    
    q0 = nf.distributions.DiagGaussian(latent_size, trainable=False)
        
    nfm = nf.NormalizingFlow(q0=q0, flows=flows)
    
    nfm = nfm.to(device)
    
    print('Model loaded')

    # --- Training --- #

    num_epochs = 100

    optimizer = torch.optim.Adam(nfm.parameters(), lr=1e-3)
    
    train_losses = []
    val_losses = []
    
    with tqdm(total=num_epochs, desc="Epochs") as pbar:
        for epoch in range(1, num_epochs + 1):
        
            nfm.train()
            train_loss = 0.0
            n_train = 0
            for x, c in train_loader:
                x = x.reshape(-1, 200*2)
                x = x.to(device)       
                c = c.to(device)
        
                loss = nfm.forward_kld(x)
        
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
        
                batch_size = x.size(0)
                train_loss += loss.item() * batch_size
                n_train += batch_size
    
            train_loss /= n_train
            train_losses.append(train_loss)
        
            nfm.eval()
            val_loss = 0.0
            n_val = 0
        
            with torch.no_grad():
                for x_val, c_val in val_loader:
                    x_val = x_val.reshape(-1, 200*2)
                    x_val = x_val.to(device)
                    c_val = c_val.to(device)
        
                    loss_val = nfm.forward_kld(x_val)
                    
                    batch_size = x_val.size(0)
                    val_loss += loss_val.item() * batch_size
                    n_val += batch_size
        
            val_loss /= n_val
            val_losses.append(val_loss)
        
            pbar.set_postfix({
                "Epoch": f"{epoch:03d}",
                "train NLL": f"{train_loss:.4f}",
                "val NLL": f"{val_loss:.4f}"
            })
            pbar.update(1)
            
    torch.save({
        "train_loss": train_losses,
        "val_loss": val_losses,
        "state_dict": nfm.state_dict(),
        "lr": optimizer.param_groups[0]['lr']
    }, f=f'/gpfs/bwfor/work/ws/hd_gy283-my_data/models/spline_flow_AR_K={K}_LS={latent_size}_HL{hidden_layers}_epochs={num_epochs}.pth'
    )