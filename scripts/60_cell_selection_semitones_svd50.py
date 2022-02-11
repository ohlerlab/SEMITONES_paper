import os
import numpy as np
import pandas as pd
import pickle
import random

from math import ceil
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.metrics.pairwise import pairwise_kernels


# load the matrix to compute similarities over
X = np.loadtxt("../data/processed/greenleaf_scRNA_combined_svd50.txt")
n_select = ceil(X.shape[0] * .002)

# select using a random starting cell
seeds = [1 * 2**i for i in range(100)]
selected = []
for seed in seeds:
    random.seed(seed)
    init = random.sample(range(X.shape[0]), 1)[0]
    selected.append(from_knn_dist(X, start=init, n_ret=n_select,
                                  metric="cosine", seed=42))

with open("cell_selection_semitones_svd50_random_init.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()

# select using the furhest cell from the medoid
S = pairwise_kernels(X, metric="rbf", gamma=.8)  # compute similarities
medoid = np.argmax(np.sum(S, axis=0))
init = np.argmin(S[:, medoid])  # furthest from the medoid
selected = from_knn_dist(X, start=init, n_ret=n_select, metric="cosine",
                         seed=42)

# save the selected cell indices
with open("cell_selection_semitones_svd50_smart_init.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()
s