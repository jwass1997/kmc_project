import torch
import torch.nn as nn

class SMPos(nn.Module):
    
    def __init__(self, num_coords, num_cvolts, layer_dims):

        super().__init__()

        self.num_coords = num_coords
        self.num_cvolts = num_cvolts

        self.layer_dims = layer_dims
        self.num_layers = len(layer_dims)
        self.in_features = layer_dims[0]
        self.out_features = layer_dims[-1]
        self.layer_list = self.build_model()

        self.model = nn.Sequential(*self.layer_list)
    
    def build_model(self):

        layers = []
        
        for k in range(1, self.num_layers):

            layers.append(nn.Linear(self.layer_dims[k-1], self.layer_dims[k]))

        return layers            


    def forward(self, x):

        out = self.model(x)

        return out