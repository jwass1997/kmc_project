import torch
import torch.nn as nn
from modules import ISAB, SAB, MAB, PMA
import torch.distributions as D

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

class CouplingNet(nn.Module):
    def __init__(self, in_dim, cond_dim, h_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim + cond_dim, h_dim),
            nn.ReLU(),
            nn.Linear(h_dim, h_dim),
            nn.ReLU()
        )
        self.s = nn.Linear(h_dim, in_dim)
        self.t = nn.Linear(h_dim, in_dim)

    def forward(self, x, c):
        h = self.net(torch.cat([x, c], dim=-1))
        s = self.s(h)
        t = self.t(h)
        
        return s, t

class CondRealNVP(nn.Module):
    def __init__(self, x_dim, cond_dim, num_coupling_layers=4, h_dim=64):
        """
        x_dim: dimension of data x (not condition)
        cond_dim: dimension of conditioning vector c
        """
        super().__init__()
        self.x_dim = x_dim
        self.cond_dim = cond_dim
        self.num_coupling_layers = num_coupling_layers

        masks = []
        # simple pattern: first half masked, then second half, alternating
        # works for any x_dim >= 2
        mask_half = torch.cat([
            #torch.ones(x_dim // 2),
            torch.ones(1),
            #torch.zeros(x_dim - x_dim // 2)
            torch.zeros(x_dim - 1)
        ])

        for i in range(num_coupling_layers):
            if i % 2 == 0:
                mask = mask_half
            else:
                mask = 1.0 - mask_half  # complement
            masks.append(mask)

        # register as buffer so it moves with .to(device) and saves in state_dict
        self.register_buffer("masks", torch.stack(masks))  # (num_layers, x_dim)

        self.coupling_layers = nn.ModuleList([
            CouplingNet(in_dim=x_dim, cond_dim=cond_dim, h_dim=h_dim)
            for _ in range(num_coupling_layers)
        ])

    def forward(self, x, c):
        """
        x: (batch, x_dim)
        c: (batch, cond_dim)
        returns: z, log_det
        """
        z = x
        log_det = torch.zeros(x.shape[0], device=x.device)

        for i, cl in enumerate(self.coupling_layers):
            mask = self.masks[i].to(z.device)      # (x_dim,)
            mask = mask.unsqueeze(0)               # (1, x_dim) for broadcasting

            x_masked = z * mask                    # part that stays
            x_trans  = z * (1.0 - mask)            # part to transform

            s, t = cl(x_masked, c)                 # (batch, x_dim)

            # only transform the unmasked dimensions
            s = s * (1.0 - mask)
            t = t * (1.0 - mask)

            z = x_masked + (x_trans * torch.exp(s) + t)

            # log-det only from transformed dims
            log_det += (s).sum(dim=-1)

        return z, log_det

    def log_prob(self, x, c):
        z, log_det = self.forward(x, c)
        base_dist = D.MultivariateNormal(
            loc=torch.zeros(self.x_dim, device=x.device),
            covariance_matrix=torch.eye(self.x_dim, device=x.device)
        )
        return base_dist.log_prob(z) + log_det

    def inverse(self, z, c):
        """
        z: (batch, x_dim)  -- sample from base_dist
        c: (batch, cond_dim)
        returns: x, log_det   (log_det is optional for sampling, but nice to keep)
        """
        x = z
        log_det = torch.zeros(z.shape[0], device=z.device)
    
        # iterate through layers in reverse
        for i in reversed(range(self.num_coupling_layers)):
            cl = self.coupling_layers[i]
            mask = self.masks[i].to(x.device).unsqueeze(0)  # (1, x_dim)
    
            x_masked = x * mask
            x_trans  = x * (1.0 - mask)
    
            # same as forward, but we condition on the masked part of *current* x
            s, t = cl(x_masked, c)
            s = s * (1.0 - mask)
            t = t * (1.0 - mask)
    
            # invert the affine transform
            x = x_masked + (x_trans - t) * torch.exp(-s)
    
            # inverse log-det is negative of forward
            log_det -= s.sum(dim=-1)
    
        return x, log_det

class CondRealNVP2d(nn.Module):
    def __init__(self, cond_dim, num_coupling_layers=4, h_dim=64):
        super().__init__()
        self.masks = []
        self.coupling_layers = nn.ModuleList()
        for i in range(num_coupling_layers):
            mask = torch.tensor([i % 2, (i + 1) % 2], dtype=torch.float32)
            self.masks.append(mask)
            self.coupling_layers.append(
                CouplingNet(in_dim=2, cond_dim=cond_dim, h_dim=h_dim)
            )
        #self.base_dist = D.MultivariateNormal(
        #    loc=torch.zeros(2).to(device),
        #    covariance_matrix=torch.eye(2).to(device)
        #)
    def forward(self, x, c):
        z = x
        log_det = torch.zeros(x.shape[0], device=x.device)
        for mask, cl in zip(self.masks, self.coupling_layers):
            mask = mask.to(z.device)
            x_masked = z * mask
            x_trans = z * (1 - mask)

            s, t = cl(x_masked, c)
            s = s * (1 - mask)
            t = t * (1 - mask)

            z = x_masked + (x_trans * torch.exp(s) + t)
            log_det += s.sum(dim=-1)
        return z, log_det
        
    def log_prob(self, x, c):
        z, log_det = self.forward(x, c)
        base_dist = D.MultivariateNormal(
            loc=torch.zeros(2, device=x.device),
            covariance_matrix=torch.eye(2, device=x.device)
        )
        return base_dist.log_prob(z) + log_det

    def inverse(self, z, c):

        x = z
        log_det = torch.zeros(x.shape[0], device=x.device)

        for i in reversed(range(self.num_coupling_layers)):
            cl = self.coupling_layers[i]
            mask = self.masks[i].to(x.device)

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