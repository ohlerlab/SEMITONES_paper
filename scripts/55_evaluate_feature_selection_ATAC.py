import os
import pickle
import numpy as np
import pandas as pd

from sklearn.neighbors import KDTree
from SEMITONES.enrichment_scoring import pvals_per_gene
from SEMITONES.tfidf import TFIDF
from SEMITONES.support_funcs import load_sparse_h5
from sklearn.decomposition import TruncatedSVD
from umap import UMAP


def density_index(X, k=None):

    """Taken from https://www.biorxiv.org/content/10.1101/2020.10.07.330563v1
    Inspired by the root-mean-squared distance drms between pairs of cells
    i and j can be calculated as."""
    
    k = 5 if k is None else k

    km = np.mean(KDTree(X).query(X, k=k)[0])
    fronorm = np.linalg.norm(X, ord="fro")
    N = X.shape[0]

    return np.sqrt(2 / N) * (fronorm / km)


def binarize(X, t):
    X = X.copy()
    condition = X.data > t
    X.data[condition] = 1
    X.data[np.logical_not(condition)] = 0
    return X


# set working directory
os.chdir("../data")
# load ATAC-matrix
X = load_sparse_h5("scATAC", "greenleaf_scATAC_peaks_filtered.h5")
X = binarize(X, 0)  # binarize the matrix

# load the peak indices
with open("greenleaf_scATAC_peaks_filtered.txt", "r") as f:
    peaks = [g.strip("\n") for g in f.readlines()]
f.close()

# compute p-value per genes
KS_per_gene = pvals_per_gene(escores, pscores).sort_values(by="q")

n_peaks = [500, 5000, 10000, 50000, 100000]  # peaks to select

DI = {}  # results dict
for n in n_peaks:
    DI[n] = {}
    
    sem_peaks = KS_per_gene.sort_values(by="q").index[:n]  # top peaks KS-test
    sem_idx = [peaks.index(i) for i in sem_peaks]  # peak indices top KS-test

    top_peaks_idx = np.argsort(-X.sum(0).A[0])[:n]  # most accessible peaks

    # load n highest ranked peaks
    filename = f"semitones_ranked_feature_selection_{n}_ATAC.pkl"
    with open(filename, "rb") as f:
        sem_rank = pickle.load(f)
    f.close()
    sem_rank_idx = [peaks.index(i) for i in sem_rank]  # get peak indices

    methods = ["semitones", "top_peaks", "sem_rank"]  # names
    indices = [sem_idx, top_peaks_idx, sem_rank_idx]  # corresponding objects
    
    for m, pset in zip(methods, indices):
        
        tfidf = TFIDF(X[:, pset])  # TF-IDF transform
        
        svd = TruncatedSVD(n_components=50).fit_transform(tfidf)  # get 50 PCs
        np.save(f"{m}_svd50_{n}peaks.npy", svd)  # save 50 PCs
        
        umap = UMAP(n_components=2, n_neighbors=50,
                    min_dist=0.5, random_state=42
                   ).fit_transform(svd)  # get UMAP
        np.save(f"{m}_umap_{n}peaks.npy", umap)  # save UMAP
        
        DI[n][m] = density_index(svd)  # compute density index

# save density indices
with open("density_indices_ATAC.pkl", "wb") as f:
    pickle.dump(DI, f)
f.close()
