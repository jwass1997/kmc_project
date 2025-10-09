import torch
import copy
import torch.nn as nn

class NeuralNet(nn.Module):

    """
    Simple Feedforward Neural Network with optional dropout on hidden layers.
    - Input/target normalization handled via registered buffers (x_mean/std, y_mean/std).
    - He (Kaiming) init for Linear layers (ReLU nonlinearity).
    """

    def __init__(
        self,
        layer_dims,
        x_mean, x_std,
        y_mean, y_std,
        dropout_p: float = 0.2,
        batch_norm: bool = False,
    ):
        super().__init__()
        self.layer_dims = layer_dims
        self.in_features = layer_dims[0]
        self.out_features = layer_dims[-1]
        self.num_layers = len(layer_dims)
        self.dropout_p = float(dropout_p)
        self.batch_norm = bool(batch_norm)

        self.register_buffer("x_mean", torch.as_tensor(x_mean, dtype=torch.float32).view(1, -1))
        self.register_buffer("x_std", torch.as_tensor(x_std, dtype=torch.float32).clamp_min(1e-8).view(1, -1))
        self.register_buffer("y_mean", torch.as_tensor(y_mean, dtype=torch.float32).view(1, -1))
        self.register_buffer("y_std", torch.as_tensor(y_std, dtype=torch.float32).clamp_min(1e-8).view(1, -1))

        self.model = self._build_model()
        self._reset_parameters()

    def _build_model(self) -> nn.Sequential:
        layers = []
        for k in range(1, len(self.layer_dims)):
            in_dim = self.layer_dims[k - 1]
            out_dim = self.layer_dims[k]

            layers.append(nn.Linear(in_dim, out_dim, bias=True))

            if k < self.num_layers - 1:
                layers.append(nn.ReLU())
                if self.dropout_p > 0:
                    layers.append(nn.Dropout(self.dropout_p))
        return nn.Sequential(*layers)

    def _reset_parameters(self):
        for m in self.model:
            if isinstance(m, nn.Linear):
                nn.init.kaiming_uniform_(m.weight, a=0.0, nonlinearity="relu")
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:

        x = (x - self.x_mean) / self.x_std
        return self.model(x)

class LearnableActivationFunc(nn.Module):
    
    def __init__(self, model, v_min, v_max, device):
        super().__init__()

        self.v_min = v_min
        self.v_max = v_max

        self.sm = copy.deepcopy(model)

        self.device = device

        if device is not None:
            self.sm = self.sm.to(device)

        for p in self.sm.parameters():
            if p.requires_grad == True:
                p.requires_grad = False
        self.sm.eval()
        
        self.reset_cntrl_parameters()

    def reset_cntrl_parameters(self):

        num_control_parameters = self.sm.layer_dims[0] - 1

        with torch.no_grad():
            init_params = (self.v_max - self.v_min) * torch.randn(num_control_parameters, device=self.device)
            #init_params = torch.tensor([1.3929118495476631, 1.4696805209710357, -1.2824866848281564, -1.4634759777883688, -0.7459804457812129, -0.2399690924094745], dtype=torch.float32).to(self.device)

        self.raw_cntrl_params = nn.Parameter(init_params)

    def scale_cntrl(self):

        return self.v_min + (self.v_max - self.v_min) * (torch.tanh(self.raw_cntrl_params) + 1) / 2

    def affine_input_scaling(self, x, eps=1e-8):

        x_min = x.min()
        x_max = x.max()

        scale = (self.v_max - self.v_min) / (x_max - x_min + eps)
        
        return self.v_min + (x - x_min) * scale 

    def forward(self, x):

        orig_shape = x.shape

        x = x.reshape(-1, 1)
        B = x.size(0)

        cntrl_params = self.scale_cntrl().to(x.device, x.dtype)[None,:].expand(B, -1)
        x_sc = self.affine_input_scaling(x)
        input_tensor = torch.cat([x_sc, cntrl_params], dim=1)

        out = self.sm(input_tensor)

        return out.reshape(orig_shape)