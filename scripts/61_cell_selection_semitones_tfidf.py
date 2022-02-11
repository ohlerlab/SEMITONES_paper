import os
import numpy as np
import pandas as pd
import pickle

from math import ceil
from scipy.io import mmread
from scipy.sparse import csr_matrix
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.tfidf import TFIDF
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.metrics.pairwise import pairwise_distances


n_cells = np.load(
    "../data/processed/greenleaf_scRNA_combined_umap35.npy"
).shape[0]

n_select = [n_cells * .001,  # .1 percent
            75,
            n_cells * .005,
            n_cells * .01]
n_select = [ceil(i) for i in n_select]

r_cells = {}

### select cells from log-normalized counts ###
print("Load TF-IDF matrix")
X = mmread("../data/interim/greenleaf_scRNA_combined_tfidf.mtx")
X = csr_matrix(X)

# using the cosine similarity
r_cells["cosine"] = {}
S = pairwise_kernels(X, metric="cosine")  # compute similarity
medoid = np.argmax(np.sum(S, axis=0))
init = np.argmin(S[:, medoid])  # furthest from the medoid
for n in n_select:
    r_cells["cosine"][n] = from_knn_dist(
        X, start=init, n_ret=n, metric="cosine", seed=42)

# using the rbf similarity
r_cells["rbf"] = {}
S = pairwise_kernels(X, metric="rbf", gamma=.8)  # compute similarity
medoid = np.argmax(np.sum(S, axis=0))
init = np.argmin(S[:, medoid])  # furthest from the medoid
for n in n_select:
    r_cells["rbf"][n] = from_knn_dist(
        X, start=init, n_ret=n, metric="rbf", metric_params={"gamma": .8},
        seed=42)

# using the euclidean distance
r_cells["euclidean"] = {}
S = pairwise_distances(X, metric="euclidean")  # compute similarity
medoid = np.argmin(np.sum(S, axis=0))
init = np.argmax(S[:, medoid])  # furthest from the medoid
for n in n_select:
    r_cells["euclidean"][n] = from_knn_dist(
        X, start=init, n_ret=n, metric="euclidean", seed=42)

# save selected cells
rname = "SEMITONES_tfidf.pkl"
with open(rname, "wb") as f:
    pickle.dump(r_cells, f)
f.close()
