import os
import numpy as np
import pandas as pd
import pickle

from math import ceil
from scipy.io import mmread
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.tfidf import TFIDF
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.metrics.pairwise import pairwise_distances


n_cells = np.loadtxt(
    "../data/processed/greenleaf_scRNA_combined_svd50.txt"
).shape[0]

# number of cells to select
n_select = [n_cells * .001,  # .1 percent
            75,
            n_cells * .005,
            n_cells * .01]
n_select = [ceil(i) for i in n_select]

### select cells from log-normalized counts ###
print("Load TF-IDF matrix")
X = mmread("../data/interim/greenleaf_scRNA_combined_tfidf.mtx")
print("Reduce to 100 dimensions")
SVD = TruncatedSVD(n_components=100).fit_transform(X)
del X

# over these dims
n_dims = [10, 30, 50, 70, 100]
for n_dim in n_dims:

    # select subset of dimensions
    X = SVD[:, :n_dim]

    # initialize result dict
    r_cells = {}
    
    # using the cosine similarity
    r_cells["cosine"] = {}
    S = pairwise_kernels(X, metric="cosine")  # compute similarity
    medoid = np.argmax(np.sum(S, axis=0))
    init = np.argmin(S[:, medoid])  # furthest from the medoid
    for n in n_select:
        r_cells["cosine"][n] = from_knn_dist(
            X, start=init, n_ret=n, metric="cosine", seed=42)

    # using the rbf kernel
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

    # save indices of selected cells
    with open(f"SEMITONES_svd{str(n_dim)}.pkl",
              "wb") as f:
        pickle.dump(r_cells, f)
    f.close()
