import torch
import torch.nn as nn
import numpy as np
import argparse

from modules import SAB, ISAB, MAB, PMA
from typing import Optional
from torch.utils.data import TensorDataset, DataLoader, random_split
from pathlib import Path
from train_sm import MakeDataset

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class PositionFFNN(nn.Module):
    
    def __init__(self, layer_dims):

        super().__init__()

        self.layer_dims = layer_dims
        self.num_layers = len(layer_dims)
        self.in_features = layer_dims[0]
        self.out_features = layer_dims[-1]

        self.model = self.build_model()
    
    def build_model(self):

        layers = []
        
        for k in range(1, self.num_layers):

            layers.append(nn.Linear(self.layer_dims[k-1], self.layer_dims[k]))

            if k < self.num_layers - 1:
                layers.append(nn.ReLU())

        return nn.Sequential(*layers)            


    def forward(self, x):

        out = self.model(x)

        return out
    
class SetTransformer(nn.Module):
    def __init__(self, dim_input, num_outputs, dim_output,
            num_inds=32, dim_hidden=128, num_heads=4, ln=False):
        super(SetTransformer, self).__init__()
        self.enc = nn.Sequential(
                ISAB(dim_input, dim_hidden, num_heads, num_inds, ln=ln),
                ISAB(dim_hidden, dim_hidden, num_heads, num_inds, ln=ln))
        self.dec = nn.Sequential(
                PMA(dim_hidden, num_heads, num_outputs, ln=ln),
                SAB(dim_hidden, dim_hidden, num_heads, ln=ln),
                SAB(dim_hidden, dim_hidden, num_heads, ln=ln),
                nn.Linear(dim_hidden, dim_output))

    def forward(self, X):
        return self.dec(self.enc(X))
    
parser = argparse.ArgumentParser()
parser.add_argument("--data_directory", type=str)
parser.add_argument("--num_batches", type=int)
parser.add_argument("--epochs", type=int)

args = parser.parse_args()

coords_list = []
voltages_list = []
I_list = []

dummy_data = np.load(Path(args.data_directory) / f"batch_{0}.npz")
input_idx = dummy_data["inputIdx"]
output_idx = dummy_data["outputIdx"]

max_voltage = 1.5
min_voltage = -1.5

for b in range(args.num_batches):
    data = np.load(Path(args.data_directory)/ f"batch_{b}.npz")
    coords_list.append(data["acc_xy"])
    voltages_list.append(data["inputs"])
    I_list.append(data["currents"])

coords_data = np.concatenate(coords_list)
voltages_data = np.concatenate(voltages_list)
voltages_data = np.delete(voltages_data, output_idx, axis=1)
I_data = np.concatenate(I_list)

N = coords_data.shape[0]
train_size = int(0.8*N)
val_size = N - train_size

coords_train = torch.from_numpy(coords_data[:N, :]).float()
voltages_train = torch.from_numpy(voltages_data[:N, :]).float()
y_train = torch.from_numpy(I_data[:N]).float()

coords_val = torch.from_numpy(coords_data[N:, :]).float()
voltages_val = torch.from_numpy(voltages_data[N:, :]).float()
y_val = torch.from_numpy(I_data[N:]).float()

y_mean = y_train.mean()
y_std = y_train.std()

y_train = (y_train - y_mean) / (y_std + 1e-12)
y_val = (y_val - y_mean) / (y_std + 1e-12)

voltages_train /= max_voltage
voltages_val /= max_voltage

train_data = TensorDataset(coords_train, voltages_train, y_train)
val_data = TensorDataset(coords_val, voltages_val, y_val)

train_loader = DataLoader(dataset=train_data, batch_size=128, shuffle=True)
val_loader = DataLoader(dataset=val_data, batch_size=128, shuffle=False)

num_points = coords_train[0].shape[0]
coords_dim = coords_train[0].shape[1]
volts_dim = voltages_train[0].shape[0]

model = SetTransformer(dim_input=coords_dim, num_outputs=1, dim_output=1)
model = model.to(device)

mse_loss = torch.nn.MSELoss()
optimzier = torch.optim.Adam(model.parameters(), lr=1e-3)

for epoch in range(1, args.epochs+1):
    model.train()
    train_loss_sum = 0.0
    for coords_, volts_, y_ in train_loader:
        coords_, volts_, y_ = coords_.to(device), volts_.to(device), y_.to(device)
        optimzier.zero_grad()
        pred = model(coords_, volts_)
        loss = mse_loss(pred, y_)
        loss.backward()
        optimzier.step()
        train_loss_sum += loss.item() * coords_.size(0)
    epoch_train_loss = train_loss_sum / len(train_loader.dataset)

    val_loss_sum = 0.0
    model.eval()
    with torch.no_grad():
        for coords_, volts_, y_ in val_loader:
            coords_, volts_, y_ = coords_.to(device), volts_.to(device), y_.to(device)
            pred = model(coords_, volts_)
            loss = mse_loss(pred, y_)
            val_loss_sum += loss.item() * coords_.size(0)
        epoch_val_loss = val_loss_sum / len(val_loader.dataset)

    print(f"[Epoch {epoch}/{args.epochs}] train_loss: {epoch_train_loss} val_loss: {epoch_val_loss}")
