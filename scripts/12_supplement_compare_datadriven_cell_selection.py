# python 3.6

import os
import gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import mmread
from SEMITONES.cell_selection import from_knn_dist
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_distances, pairwise_kernels
from SEMITONES.tfidf import TFIDF

os.chdir("../data/processed/")

### 25D UMAP ###
UMAP25 = np.load("greenleaf_scRNA_combined_umap25.npy")

# RBF over UMAP
i = np.argmin(np.sum(pairwise_kernels(UMAP25, metric="rbf",
                                      **{"gamma": 0.08}), axis=1))
UMAP_RBF = from_knn_dist(UMAP25, start=i, n_ret=75,
                         metric="rbf", seed=42,
                         metric_params={"gamma": 0.8})

# Cosine over UMAP
i = np.argmin(np.sum(pairwise_kernels(UMAP25, metric="cosine"), axis=1))
UMAP_cosine = from_knn_dist(UMAP25, start=i, n_ret=75,
                            metric="cosine", seed=42)

# Euclidean over UMAP
i = np.argmax(np.sum(pairwise_distances(UMAP25, metric="euclidean"), axis=1))
UMAP_euclidean = from_knn_dist(UMAP25, start=i, n_ret=75,
                               metric="euclidean", seed=42)

del UMAP25
gc.collect()

### 50 PCs ###
SVD = mmread("greenleaf_scRNA_combined_svd50.mtx")

# RBF over SVD
i = np.argmin(np.sum(pairwise_kernels(SVD, metric="rbf",
                                      **{"gamma": 8e-1}), axis=1))
SVD_RBF = from_knn_dist(SVD, start=i, n_ret=75,
                        metric="rbf", metric_params={"gamma": 8e-1}, seed=42)

# Cosine over SVD
i = np.argmin(np.sum(pairwise_kernels(SVD, metric="cosine"), axis=1))
SVD_cosine = from_knn_dist(SVD, start=i, n_ret=75, metric="cosine", seed=42)

# Euclidean over SVD
i = np.argmax(np.sum(pairwise_distances(SVD, metric="euclidean"), axis=1))
SVD_euclidean = from_knn_dist(SVD, start=i, n_ret=75, metric="euclidean", seed=42)

del SVD
gc.collect()

### NORM ###
NORM = mmread("greenleaf_scRNA_combined_norm.mtx")
NORM = NORM.tocsr()

# RBF over SVD
i = np.argmin(np.sum(pairwise_kernels(NORM, metric="rbf",
                                      **{"gamma": 8e-1}), axis=1))
NORM_RBF = from_knn_dist(NORM, start=i, n_ret=75,
                         metric="rbf", metric_params={"gamma": 8e-1}, seed=42)


# Cosine over SVD
i = np.argmin(np.sum(pairwise_kernels(NORM, metric="cosine"), axis=1))
NORM_cosine = from_knn_dist(NORM, start=i, n_ret=75, metric="cosine", seed=42)

# Euclidean over SVD
i = np.argmax(np.sum(pairwise_distances(NORM, metric="euclidean"), axis=1))
NORM_euclidean = from_knn_dist(NORM, start=i, n_ret=75, metric="euclidean", seed=42)

### TFIDF ###
tfidf = TFIDF(NORM)
del NORM
gc.collect()

# RBF over SVD
i = np.argmin(np.sum(pairwise_kernels(tfidf, metric="rbf",
                                      **{"gamma": 8e-1}), axis=1))
tfidf_RBF = from_knn_dist(tfidf, start=i, n_ret=75,
                          metric="rbf", metric_params={"gamma": 8e-1}, seed=42)


# Cosine over SVD
i = np.argmin(np.sum(pairwise_kernels(tfidf, metric="cosine"), axis=1))
tfidf_cosine = from_knn_dist(tfidf, start=i, n_ret=75, metric="cosine", seed=42)

# Euclidean over SVD
i = np.argmax(np.sum(pairwise_distances(tfidf, metric="euclidean"), axis=1))
tfidf_euclidean = from_knn_dist(tfidf, start=i, n_ret=75, metric="euclidean", seed=42)

del tfidf
gc.collect()

selections = [UMAP_RBF, UMAP_cosine, UMAP_euclidean,
              SVD_RBF, SVD_cosine, SVD_euclidean,
              NORM_RBF, NORM_cosine, NORM_euclidean,
              tfidf_RBF, tfidf_cosine, tfidf_euclidean]
index = ["UMAP_RBF", "UMAP_cosine", "UMAP_euclidean",
         "SVD_RBF", "SVD_cosine", "SVD_euclidean",
         "NORM_RBF", "NORM_cosine", "NORM_euclidean",
         "tfidf_RBF", "tfidf_cosine", "tfidf_euclidean"]
selections = pd.DataFrame(selections, index=index)
selections.to_csv("supplemental_cell_selection.txt", sep="\t")
