import torch
import torch.nn as nn

""" Simple Feedforward Neural Network """

class NeuralNet(nn.Module):
    
    def __init__(
            self,
            layer_dims,
            x_mean, x_std,
            y_mean, y_std,
            dropout_p: float = 0.2,
            batch_norm: bool = False
    ):
        super(NeuralNet, self).__init__()

        self.layer_dims = layer_dims
        self.in_features = layer_dims[0]
        self.out_features = layer_dims[-1]
        self.num_layers = len(layer_dims)
        self.dropout_p    = dropout_p
        self.batch_norm = batch_norm

        self.register_buffer('x_mean', torch.as_tensor(x_mean, dtype=torch.float32))
        self.register_buffer('x_std', torch.as_tensor(x_std, dtype=torch.float32))

        self.register_buffer('y_mean', torch.as_tensor(y_mean, dtype=torch.float32))
        self.register_buffer('y_std', torch.as_tensor(y_std, dtype=torch.float32))

        self.model = self.build_model()
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:

        x = (x - self.x_mean) / self.x_std
        
        out = self.model(x)

        return out
    
    def build_model(self):
        layer_list = []

        for k in range(1, len(self.layer_dims)):
            layer_list.append(nn.Linear(self.layer_dims[k-1], self.layer_dims[k]))

            if k < self.num_layers - 1:
                #if self.batch_norm:
                #    layer_list.append(nn.BatchNorm1d(self.layer_dims[k]))
                #layer_list.append(nn.ReLU(inplace=True))
                layer_list.append(nn.PReLU())

        """for l in range(self.num_layers):
            # Determine dims for this layer
            if l == 0:
                in_dim, out_dim = self.in_features, self.hidden_dim
            elif l == self.num_layers - 1:
                in_dim, out_dim = self.hidden_dim, self.out_features
            else:
                in_dim = out_dim = self.hidden_dim

            # Linear
            layer_list.append(nn.Linear(in_dim, out_dim))

            # If not final layer, add BN, activation, dropout
            if l < self.num_layers - 1:
                layer_list.append(nn.BatchNorm1d(out_dim))
                layer_list.append(nn.ReLU())
                #layer_list.append(nn.Dropout(self.dropout_p))"""

        return nn.Sequential(*layer_list)