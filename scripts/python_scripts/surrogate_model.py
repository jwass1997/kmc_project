import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

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
    
class NeuralNet(nn.Module):
    
    def __init__(
            self,
            in_features: int,
            out_features: int,
            hidden_dim: int,
            num_layers: int
    ):
        super(NeuralNet, self).__init__()

        self.in_features = in_features
        self.out_features = out_features
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers

        self.model_layers = self.build_model()
        self.model = nn.Sequential(*self.model_layers)
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        
        out = self.model(x)

        return out
    
    def build_model(self):

        layer_list = []

        for l in range(self.num_layers):

            if l == 0:
                layer_list += [nn.Linear(in_features=self.in_features, out_features=self.hidden_dim)]
                layer_list += [nn.ReLU()]

            if l == self.num_layers - 1:
                layer_list += [nn.Linear(in_features=self.hidden_dim, out_features=self.out_features)]

            else:
                layer_list += [nn.Linear(in_features=self.hidden_dim, out_features=self.hidden_dim)]
                layer_list += [nn.ReLU()]
        
        return layer_list
    
if __name__ == "__main__":

    data_0 = np.load(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/test_batch/batch_steps=1e6_0.npz")
    data_1 = np.load(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/test_batch/batch_steps=1e6_1.npz")
    data_2 = np.load(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/test_batch/batch_steps=1e6_2.npz")
    data_3 = np.load(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/test_batch/batch_steps=1e6_3.npz")
    data_4 = np.load(f"/gpfs/bwfor/work/ws/hd_gy283-my_data/test_batch/batch_steps=1e6_4.npz")

    inputs_0 = data_0["inputs"]
    inputs_1 = data_1["inputs"]
    inputs_2 = data_2["inputs"]
    inputs_3 = data_3["inputs"]
    inputs_4 = data_4["inputs"]

    outputs_0 = data_0["currents"]
    outputs_1 = data_1["currents"]
    outputs_2 = data_2["currents"]
    outputs_3 = data_3["currents"]
    outputs_4 = data_4["currents"]


    inputs = np.concatenate([inputs_0, inputs_1, inputs_2, inputs_3, inputs_4])
    outputs = np.concatenate([outputs_0, outputs_1, outputs_2, outputs_3, outputs_4])

    print(f"input_shape = {inputs.shape}")
    print(f"output_shape = {outputs.shape}")

    X_train, X_test, y_train, y_test = train_test_split(inputs, outputs, test_size=0.2, random_state=42, shuffle=True)

    train_set = MakeDataset(X_train, y_train)
    test_set = MakeDataset(X_test, y_test)

    batch_size = 128

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True,  num_workers=4)
    test_loader  = DataLoader(test_set,  batch_size=batch_size, shuffle=False, num_workers=4)

    learning_rate = 1e-5
    num_epochs = 3000

    model = NeuralNet(in_features=8, out_features=1, hidden_dim=90, num_layers=5).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    for epoch in range(1, num_epochs+1):

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

        print(f"Epoch {epoch:2d}/{num_epochs}   "
            f"Train Loss: {epoch_train_loss:.10f}   "
            f"Val Loss: {epoch_val_loss:.10f}")
    
    torch.save(model.state_dict(), "SM_num_layers=5_hd=90_0.pth")
