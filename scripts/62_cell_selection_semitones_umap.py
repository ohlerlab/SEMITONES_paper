import os
import numpy as np
import pandas as pd
import pickle

from math import ceil
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.tfidf import TFIDF
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.metrics.pairwise import pairwise_distances
from umap import UMAP


n_cells = np.loadtxt(
    "../data/processed/greenleaf_scRNA_combined_umap25.txt"
).shape[0]

n_select = [n_cells * .001,  # .1 percent
            75,
            n_cells * .005,
            n_cells * .01]
n_select = [ceil(i) for i in n_select]

### select cells from nD UMAPs ###
SVD = np.loadtxt("../data/processed/greenleaf_scRNA_combined_svd50.txt")
n_dims = [5, 10, 25, 35, 50]
for n_dim in n_dims:

    if n_dim == 25:  # take precomputed for speed
        X = np.loadtxt("../data/processed/greenleaf_scRNA_combined_umap25.txt")
    else:
        X = UMAP(n_components=n_dim, min_dist=.3, n_neighbors=30
                ).fit_transform(SVD)

    r_cells = {}  # initialize dict

    # select cells using the cosine similarity
    r_cells["cosine"] = {}
    S = pairwise_kernels(X, metric="cosine")  # compute similarities
    medoid = np.argmax(np.sum(S, axis=0))
    init = np.argmin(S[:, medoid])  # furthest from medoid
    for n in n_select:
        r_cells["cosine"][n] = from_knn_dist(
            X, start=init, n_ret=n, metric="cosine", seed=42)

    # using the rbf similarity
    r_cells["rbf"] = {}
    S = pairwise_kernels(X, metric="rbf", gamma=.8)  # compute similarities
    medoid = np.argmax(np.sum(S, axis=0))
    init = np.argmin(S[:, medoid])  # furthest from medoid
    for n in n_select:
        r_cells["rbf"][n] = from_knn_dist(
            X, start=init, n_ret=n, metric="rbf", metric_params={"gamma": .8},
            seed=42)

    # using the euclidean distance
    r_cells["euclidean"] = {}
    S = pairwise_distances(X, metric="euclidean")  # compute similarities
    medoid = np.argmin(np.sum(S, axis=0))
    init = np.argmax(S[:, medoid])  # furthest from medoid
    for n in n_select:
        r_cells["euclidean"][n] = from_knn_dist(
            X, start=init, n_ret=n, metric="euclidean", seed=42)

    # save selected cells dict
    rname = f"SEMITONES_umap{str(n_dim)}.pkl"
    with open(rname, "wb") as f:
        pickle.dump(r_cells, f)
    f.close()
