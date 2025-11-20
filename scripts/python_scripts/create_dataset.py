import torch
import torch.nn as nn
import networkx as nx
import json
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split
from train_sm import train_val_test_loaders, MakeDataset
from pathlib import Path

data_dir = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/data_uni_samples_1e6/')
save_dir = Path('/gpfs/bwfor/work/ws/hd_gy283-my_data/data/')
save_dir.mkdir(exist_ok=True, parents=True)

n_a = 200
radius = 150.0
R = np.sqrt(np.pi * radius**2 / n_a)

output_idx = 0
num_batches = 200

if __name__ == '__main__':

    for i in range(num_batches):

        fp = data_dir / f'batch_{i}.npz'
        try:
            batch = np.load(fp)
        except FileNotFoundError:
            print(f"Batch file {fp} not found, skipping.")
            continue

        # --- extract data for this batch only ---
        V = batch['inputs']       # shape: (B, num_inputs)
        X = batch['acc_xy']       # shape: (B, N_a, 2)
        A_list = [a_ij[:n_a, :n_a] for a_ij in batch['adj_mat']]

        if batch['currents'].ndim != 1:
            y = batch['currents'][:, output_idx]
        else:
            y = batch['currents']

        # drop output_idx column from V
        keep_cols = [j for j in range(V.shape[1]) if j != output_idx]
        V = V[:, keep_cols]

        # rescale positions
        X = X * R / radius

        # stack adjacency matrices into a single array for this batch
        A = np.stack(A_list, axis=0)

        # optionally downcast to save memory & disk (if okay for you):
        V = V.astype(np.float32)
        X = X.astype(np.float32)
        A = A.astype(np.float32)       # or np.uint8 if adjacency is 0/1
        y = y.astype(np.float32)

        # --- compute graph stats per sample for this batch ---
        num_samples = A.shape[0]
        avg_degrees = np.empty(num_samples, dtype=np.float32)
        avg_clustering = np.empty(num_samples, dtype=np.float32)

        for k, a_ij in enumerate(A):
            num_nodes = a_ij.shape[0]
            num_edges = np.count_nonzero(np.triu(a_ij, k=1))

            avg_deg_k = 2.0 * num_edges / num_nodes
            avg_degrees[k] = avg_deg_k

            g = nx.from_numpy_array(a_ij)
            avg_clst_k = nx.average_clustering(g)
            avg_clustering[k] = avg_clst_k

            if k % 1000 == 0:
                print(f"Batch {i}, sample {k}/{num_samples}")

        # --- save this batch as its own dataset file ---
        out_path = save_dir / f'dataset_batch_{i:04d}.pt'
        torch.save(
            {
                "V": V,
                "X": X,
                "A": A,
                "avg_degrees": avg_degrees,
                "avg_clustering": avg_clustering,
                "y": y,
                "radius": radius,
                "out_idx": output_idx,
            },
            out_path
        )
        print(f"Saved {out_path}")