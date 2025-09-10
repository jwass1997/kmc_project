import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
import argparse
from pathlib import Path
from argparse import Namespace
from sm import NeuralNet

from sklearn.model_selection import train_test_split

class MakeDataset(Dataset):
    def __init__(self, x_data, labels=None):

        if isinstance(x_data, torch.Tensor):
            self.X = x_data.float()
        elif isinstance(x_data, np.ndarray):
            self.X = torch.from_numpy(x_data).float()
        else:
            raise TypeError(f"x_data must be a numpy array or torch.Tensor, got {type(x_data)}")

        if labels is not None:
            if isinstance(labels, torch.Tensor):
                y = labels.float()
            elif isinstance(labels, np.ndarray):
                y = torch.from_numpy(labels).float()
            else:
                raise TypeError(f"labels must be a numpy array or torch.Tensor, got {type(labels)}")

            if y.dim() == 1:
                y = y.unsqueeze(1)

            if y.size(0) != self.X.size(0):
                raise ValueError("inputs and labels must have the same length")

            self.y = y
        else:
            self.y = None

    def __len__(self):
        return self.X.size(0)

    def __getitem__(self, idx):
        x = self.X[idx]
        if self.y is None:
            return x
        return x, self.y[idx]
    
def create_data_loaders(data_dir, num_batches, batch_size, normalize=False, eps=1e-8, train_size=0.8, test_size=0.2, random_state=42, num_workers=0):

    input_list = []
    output_list = []

    for i in range(num_batches):
        batch = np.load(data_dir/f"batch_1e6_{i}.npz")
        _input = batch["inputs"]
        _output = batch["currents"]
        input_list.append(_input)
        output_list.append(_output)

    raw_inputs = np.concatenate(input_list)
    outputs = np.concatenate(output_list)
    inputs = raw_inputs[:, 1:]

    print(f"raw_inputs shape: {raw_inputs.shape}")
    print(f"outputs shape: {outputs.shape}")
    print(f"inputs shape: {inputs.shape}")

    X_train, X_test, y_train, y_test = train_test_split(
        inputs, 
        outputs, 
        test_size=test_size, 
        random_state=random_state, 
        shuffle=True
    )

    if normalize:
        X_train_mean = X_train.mean(0)
        X_train_std = X_train.std(0) + eps

        y_train_mean = y_train.mean(0)
        y_train_std = y_train.std(0) + eps

        X_train = (X_train - X_train_mean) / X_train_std
        X_test = (X_test - X_train_mean) / X_train_std
        y_train = (y_train - y_train_mean) / y_train_std
        y_test= (y_test - y_train_mean) / y_train_std
    
    train_set = MakeDataset(X_train, y_train)
    test_set = MakeDataset(X_test, y_test)

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True, num_workers=num_workers, pin_memory=True)
    test_loader  = DataLoader(test_set,  batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)

    return train_loader, test_loader      

def train_sm(model, criterion, optimizer, args, device, train_loader, test_loader, print_every=100):

    train_losses = []
    val_losses = []

    for epoch in range(1, args.num_epochs+1):

        model.train()
        running_loss = 0.0
        for inputs_batch, targets_batch in train_loader:
            inputs_batch  = inputs_batch.to(device, non_blocking=True)
            targets_batch = targets_batch.to(device, non_blocking=True)

            optimizer.zero_grad()
            preds = model(inputs_batch)
            loss  = criterion(preds, targets_batch)
            loss.backward()
            optimizer.step()

            running_loss += loss.item() * inputs_batch.size(0)

        epoch_train_loss = running_loss / len(train_loader.dataset)

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for inputs_batch, targets_batch in test_loader:
                inputs_batch  = inputs_batch.to(device, non_blocking=True)
                targets_batch = targets_batch.to(device, non_blocking=True)

                preds = model(inputs_batch)
                loss  = criterion(preds, targets_batch)
                val_loss += loss.item() * inputs_batch.size(0)

        epoch_val_loss = val_loss / len(test_loader.dataset)

        train_losses.append(epoch_train_loss)
        val_losses.append(epoch_val_loss)

        if epoch % print_every == 0:
            print(f"Epoch {epoch:2d}/{args.num_epochs} "
                f"   Train Loss: {epoch_train_loss:.10f}"
                f"   Val Loss: {epoch_val_loss:.10f}", flush=True)

    torch.save({
        "model_state_dict": model.state_dict(),
        "args": vars(args),
        "train_losses": train_losses,
        "val_losses": val_losses
    }, f"{args.save_name}.pth")

if __name__ == "__main__":

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    num_batch_list = [10, 20, 50, 100]

    criterion = nn.MSELoss()

    for num_batches in num_batch_list:

        args = Namespace(
            data_dir=Path(f"/home/hd/hd_hd/hd_gy283/kmc_project/data/sm_batches_1e6_vonMises_beta"),
            num_batches=num_batches,
            batch_size=128,
            normalize=False,
            hd=90,
            num_layers=8,
            num_epochs=1_000,
            lr=1e-3,
            save_name=Path(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/sm_vMB_num_batches={num_batches}"),
        )

        train_loader, test_loader = create_data_loaders(args.data_dir, args.num_batches, args.batch_size, args.normalize)

        model = NeuralNet(in_features=7, out_features=1, hidden_dim=args.hd, num_layers=args.num_layers)
        model = model.to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

        train_sm(model=model, criterion=criterion, optimizer=optimizer, args=args, device=device, train_loader=train_loader, test_loader=test_loader)