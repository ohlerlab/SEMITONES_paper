import numpy as np
import os
import pickle
import random

from math import ceil


# laod the 50D SVD
svd_file = f"../data/processed/greenleaf_scRNA_combined_svd50.txt"
SVD = np.loadtxt(svd_file)

geom_prog = [1 * 2**i for i in range(100)]  # seeds

n = ceil(SVD.shape[0] * 0.002)  # .2% of cells
selected = []
for seed in geom_prog:
    random.seed(seed)
    selected.append(random.sample(range(SVD.shape[0]), n))

# save selected cell indices
with open("cell_selection_random_sampling.pkl",
          "wb") as f:
    pickle.dump(selected, f)
f.close()
