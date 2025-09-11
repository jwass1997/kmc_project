import torch
import torch.nn as nn

""" Simple Feedforward Neural Network """

class NeuralNet(nn.Module):
    
    def __init__(
            self,
            in_features: int,
            out_features: int,
            hidden_dim: int,
            num_layers: int,
            dropout_p: float = 0.2
    ):
        super(NeuralNet, self).__init__()

        self.in_features = in_features
        self.out_features = out_features
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout_p    = dropout_p

        self.model_layers = self.build_model()
        self.model = nn.Sequential(*self.model_layers)
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        
        out = self.model(x)

        return out
    
    def build_model(self):
        layer_list = []

        for l in range(self.num_layers):
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
                #layer_list.append(nn.Dropout(self.dropout_p))

        return layer_list