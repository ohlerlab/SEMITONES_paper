import geosketch
import numpy as np
import os
import pickle

from math import ceil


# embedding to select cells from
svd_file = "../data/processed/greenleaf_scRNA_combined_svd50.txt"
SVD = np.loadtxt(svd_file)

n = ceil(SVD.shape[0] * 0.002)  # select .2% of cells

# select cells using standard params
seeds = np.linspace(1, 1e6, 100, dtype="int32")
selected = []
for seed in seeds:
    selected.append(geosketch.gs(SVD, N=n, seed=seed))

# save resulting cell indices
with open("cell_selection_geosketch_default.pkl", "wb") as f:
    pickle.dump(selected, f)
f.close()

# select cells with higher k for better performance
seeds = np.linspace(1, 1e6, 100, dtype="int32")
selected = []
for seed in seeds:
    selected.append(geosketch.gs(SVD, N=n, seed=seed, k=10000))

# save resulting cell indices
with open("cell_selection_geosketch_k10000.pkl", "wb") as f:
    pickle.dump(selected, f)
f.close()
