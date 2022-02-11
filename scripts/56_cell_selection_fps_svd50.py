import numpy as np
import os
import pickle
import random

from math import ceil
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels


def greedy_fps(D, n, smart_init=None, seed=None):
    
    """Greedy farthest point sampling for reference cell selection.
    -----
    D: matrix-like
        A symmetrical matrix with pairwise distances between cells.
    n: int
        The number of reference cells to select
    smart_init: bool
        If True, select the furthest cell from the medoid as the
        starting cell. Otherwise pick a random cell.
    -----
    Returns a list of reference cell indices.
    """
    
    smart_init = True if smart_init is None else smart_init
    
    selected = []

    if smart_init:
        medoid = np.argmin(np.sum(D, axis=1))
        # initialize with furthest point from the medoid
        selected.append(np.argmax(D[:, medoid]))
    else:
        random.seed(seed)
        selected.append(random.sample(range(D.shape[0]), 1)[0])
    
    # iteratively select the next furthest cell
    selected.append(np.argmax(D[:, selected]))
    for i in range(1, n):
        selected.append(np.argmax(np.mean(D[:, selected], axis=1)))

    return selected


# select from SVD space
X = np.loadtxt("../data/processed/greenleaf_scRNA_combined_svd50.txt")
n_cells = X.shape[0]  # number of cells
n = ceil(n_cells * 0.002)  # .2% of cells

# select with random starting cell
geom_prog = [1 * 2**i for i in range(100)]  # seeds

selected = []
D = 1 - pairwise_kernels(X, metric="cosine")  # distance
for seed in geom_prog:
    selected.append(greedy_fps(D, n, smart_init=False, seed=seed))

with open("cell_selection_fps_random_init_svd50.pkl", "wb") as f:
    pickle.dump(selected, f)
f.close()

# with selecting the further point from the medoid as starting cell
selected = greedy_fps(D, n, smart_init=True)
with open("cell_selection_fps_smart_init_svd50.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()
