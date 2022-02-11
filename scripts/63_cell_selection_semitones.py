import os
import numpy as np
import pandas as pd
import pickle
import random

from math import ceil
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.metrics.pairwise import pairwise_kernels


# get the UMAP to compute over
X = np.loadtxt("../data/processed/greenleaf_scRNA_combined_umap25.txt")
n_select = ceil(X.shape[0] * .002)  # number of cells to select

seeds = [1 * 2**i for i in range(100)]

# do cell seelction with random starting cells
selected = []
for seed in seeds:
    random.seed(seed)
    init = random.sample(range(X.shape[0]), 1)[0]
    selected.append(from_knn_dist(X, start=init, n_ret=n_select, metric="rbf",
                                  metric_params={"gamma": .8}, seed=42))

# save randomly selected indices
with open("cell_selection_semitones_umap25_random_init.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()

# select furthest cell from medoid as reference cell
S = pairwise_kernels(X, metric="rbf", gamma=.8)
medoid = np.argmax(np.sum(S, axis=0))
init = np.argmin(S[:, medoid])
selected = from_knn_dist(X, start=init, n_ret=n_select, metric="rbf",
                         metric_params={"gamma": .8}, seed=42)

# save cell indices for furthest cell from medoid
with open("cell_selection_semitones_umap25_smart_init.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()
