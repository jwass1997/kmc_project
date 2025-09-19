import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from typing import Optional

class MAB(nn.Module):
    def __init__(self, dim_Q, dim_KV, dim_out, num_heads=4, ln=True):
        super().__init__()
        self.num_heads = num_heads
        self.fc_q = nn.Linear(dim_Q, dim_out)
        self.fc_k = nn.Linear(dim_KV, dim_out)
        self.fc_v = nn.Linear(dim_KV, dim_out)
        self.fc_o = nn.Linear(dim_out, dim_out)
        self.ln0 = nn.LayerNorm(dim_out) if ln else nn.Identity()
        self.ln1 = nn.LayerNorm(dim_out) if ln else nn.Identity()
        self.ff = nn.Sequential(
            nn.Linear(dim_out, 4*dim_out), nn.GELU(), nn.Linear(4*dim_out, dim_out)
        )


def forward(self, Q, K, V, mask_Q: Optional[torch.Tensor]=None, mask_K: Optional[torch.Tensor]=None):
    # Q: (B,nQ,dq), K,V: (B,nKV,dkv)
    Qh, Kh, Vh = self.fc_q(Q), self.fc_k(K), self.fc_v(V)
    B, nQ, d = Qh.shape
    _, nKV, _ = Kh.shape
    h = self.num_heads
    d_h = d // h
    def split(x, n):
        return x.view(B, n, h, d_h).transpose(1, 2) # (B,h,n,d_h)
    Qm, Km, Vm = split(Qh, nQ), split(Kh, nKV), split(Vh, nKV)
    att = (Qm @ Km.transpose(-2, -1)) / math.sqrt(d_h) # (B,h,nQ,nKV)
    if mask_K is not None:
        mask = mask_K[:, None, None, :].expand(B, h, nQ, nKV)
        att = att.masked_fill(~mask, float('-inf'))
    A = att.softmax(dim=-1)
    H = A @ Vm # (B,h,nQ,d_h)
    H = H.transpose(1, 2).contiguous().view(B, nQ, d)
    H = self.ln0(Qh + self.fc_o(H))
    H = self.ln1(H + self.ff(H))
    return H


class SAB(nn.Module):
    def __init__(self, dim_in, dim_out, num_heads=4, ln=True):
        super().__init__()
        self.mab = MAB(dim_in, dim_in, dim_out, num_heads, ln)
    def forward(self, X, mask=None):
        return self.mab(X, X, X, mask_Q=mask, mask_K=mask)


class ISAB(nn.Module):
    def __init__(self, dim_in, dim_out, m=32, num_heads=4, ln=True):
        super().__init__()
        self.I = nn.Parameter(torch.randn(1, m, dim_out))
        self.mab1 = MAB(dim_out, dim_in, dim_out, num_heads, ln)
        self.mab2 = MAB(dim_in, dim_out, dim_out, num_heads, ln)
    def forward(self, X, mask=None):
        B = X.size(0)
        H = self.mab1(self.I.expand(B, -1, -1), X, X, mask_K=mask)
        return self.mab2(X, H, H)


class PMA(nn.Module):
    def __init__(self, dim_in, k=1, num_heads=4, ln=True):
        super().__init__()
        self.S = nn.Parameter(torch.randn(1, k, dim_in))
        self.mab = MAB(dim_in, dim_in, dim_in, num_heads, ln)
    def forward(self, X, mask=None):
        B = X.size(0)
        return self.mab(self.S.expand(B, -1, -1), X, X, mask_K=mask)