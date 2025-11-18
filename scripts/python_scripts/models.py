import torch
import torch.nn as nn
from modules import ISAB, SAB, MAB, PMA

class DeepSet(nn.Module):
    def __init__(self, in_dim, h_dim, num_outputs, out_dim, aggr_type):
        super().__init__()

        self.in_dim = in_dim
        self.h_dim = h_dim
        self.num_outputs = num_outputs
        self.out_dim = out_dim

        self.rho = nn.Sequential(
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, out_dim*num_outputs),
        )
        self.phi = nn.Sequential(
            nn.Linear(in_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim)
        )

        self.aggr_type = aggr_type

    def forward(self, X, mask=None):

        h = self.phi(X)
        
        if mask is not None:
            h = h * mask.unsqueeze(-1)

        if self.aggr_type == 'sum':
            s = h.sum(dim=-2)
        elif self.aggr_type == 'max':
            s, _ = h.max(dim=-2)
        elif self.aggr_type == 'mean':
            s = h.mean(dim=-2)
        else:
            raise ValueError('Unknown aggregation type')

        return self.rho(s).reshape(-1, self.num_outputs, self.out_dim)

class SetTransformer(nn.Module):
    def __init__(self, in_dim, num_outputs, out_dim,
            num_inds=32, h_dim=128, num_heads=4, ln=False):
        super(SetTransformer, self).__init__()
        self.enc = nn.Sequential(
                ISAB(in_dim, h_dim, num_heads, num_inds, ln=ln),
                ISAB(h_dim, h_dim, num_heads, num_inds, ln=ln))
        self.dec = nn.Sequential(
                PMA(h_dim, num_heads, num_outputs, ln=ln),
                SAB(h_dim, h_dim, num_heads, ln=ln),
                SAB(h_dim, h_dim, num_heads, ln=ln),
                nn.Linear(h_dim, out_dim))

    def forward(self, X):
        return self.dec(self.enc(X))

class CondSM(nn.Module):

    def __init__(
        self,
        layer_dims,
        dropout_p: float = 0.2,
        batch_norm: bool = True,
    ):
        super().__init__()

        self.layer_dims = layer_dims
        self.in_features = layer_dims[0]
        self.out_features = layer_dims[-1]
        self.num_layers = len(layer_dims)

        self.dropout_p = float(dropout_p)
        self.batch_norm = bool(batch_norm)

        self.model = self._build_model()
        self._reset_parameters()

    def _build_model(self) -> nn.Sequential:
        layers = []
        for k in range(1, len(self.layer_dims)):
            in_dim = self.layer_dims[k - 1]
            out_dim = self.layer_dims[k]

            layers.append(nn.Linear(in_dim, out_dim, bias=True))

            if k < self.num_layers - 1:
                if self.dropout_p > 0:
                    layers.append(nn.Dropout(self.dropout_p))
                if self.batch_norm:
                    layers.append(nn.BatchNorm1d(out_dim))
                layers.append(nn.ReLU())
                
        return nn.Sequential(*layers)

    def _reset_parameters(self):
        for m in self.model:
            if isinstance(m, nn.Linear):
                nn.init.kaiming_uniform_(m.weight, a=0.0, nonlinearity="relu")
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:

        return self.model(x)

"""class DeepSet(nn.Module):
    def __init__(self, in_dim, h_dim, out_dim, aggr_type):
        super().__init__()

        self.phi = nn.Sequential(
            nn.Linear(in_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim)
        )
        self.rho = nn.Sequential(
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, out_dim)
        )

        self.aggr_type = aggr_type

    def forward(self, x, mask=None):

        h = self.phi(x)

        if mask is not None:
            h = h * mask.unsqueeze(-1)
            
        if self.aggr_type == 'sum':
            s = h.sum(dim=1)
        elif self.aggr_type == 'max':
            s, _ = h.max(dim=1)
        elif self.aggr_type == 'mean':
            s = h.mean(dim=1)
        else:
            raise ValueError('Unknown aggregation type')

        return self.rho(s)"""

"""class SetTransformer(nn.Module):
    def __init__(
        self,
        dim_input=2,
        num_outputs=1,
        dim_output=6,
        num_inds=16,
        dim_hidden=128,
        num_heads=4,
        ln=False,
    ):
        super(SetTransformer, self).__init__()
        self.enc = nn.Sequential(
            ISAB(dim_input, dim_hidden, num_heads, num_inds, ln=ln),
            ISAB(dim_hidden, dim_hidden, num_heads, num_inds, ln=ln),
        )
        self.dec = nn.Sequential(
            nn.Dropout(),
            PMA(dim_hidden, num_heads, num_outputs, ln=ln),
            nn.Dropout(),
            nn.Linear(dim_hidden, dim_output),
        )

    def forward(self, X):
        return self.dec(self.enc(X)).squeeze()"""