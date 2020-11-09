# miniconda3/envs/hcv1

import os
import gc
import pickle
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.enrichment_scoring import calculate_escores, permute
gc.enable()  # enable garbage collection
os.chdir("../data/processed/")  # set data directory

def _binarize(X, t=None):
    t = 0 if t is None else t
    X = X.copy()
    condition = X.data > t
    X.data[condition] = 1
    X.data[np.logical_not(condition)] = 0
    return X


### LOAD CELL X PEAK ###
print("Loading peak accesibility data")
fname = "greenleaf_scATAC_peaks_rare.mtx"
peaks = mmread(fname)  # read in the sparse matrix
peaks = peaks.tocsr()


### BINARIZE THE DATA ###
print("Binarizing peaks")
peaks = _binarize(peaks, t=0)


### REFERENCE CELLS ###
print("Obtaining reference cells")
fname = "greenleaf_scATAC_cells_selected_from_knn.txt"
with open(fname, "r") as f:
    r_cells = [int(c.strip("\n")) for c in f.readlines()]
f.close()


### SIMILARITY TO R_CELLS ###
print("Calculating similarities to the reference cells")
# load 35D UMAP
umap35 = np.load("../interim/greenleaf_scATAC_filtered_umap35.npy")
metric = "rbf"  # radial basis function kernel
gamma = 8e-1  # influence radius
S = pairwise_similarities(umap35, r_cells, metric="rbf",
                          metric_params={"gamma": gamma})
del umap35
gc.collect()


### ENRICHMENT SCORING ###
print("Start scoring")
scores = calculate_escores(peaks, query=r_cells, S=S,
                           parallelize=True, scale_exp=True)
fname = "greenleaf_scATAC_peaks_rare_labels.txt"
with open(fname, "r") as f:
    peaknames = [p.strip("\n") for p in f.readlines()]
f.close()
scores.index = peaknames
fname = "greenleaf_scATAC_escores_knn_rbf_8e-1_rare.txt"
scores.to_csv(fname, sep="\t")
print("Saved scores")


### PERMUTE DATAFRAME ###
print("Permuting dataframe")
P = permute(peaks, n=256, axis=0, seed=42)
del peaks
gc.collect()


### P ENRICHMENT SCORING ###
print("Enrichment scoring over permutated dataframe")
pscores = calculate_escores(P, query=r_cells, S=S,
                            parallelize=True, scale_exp=True)
fname = "../interim/greenleaf_scATAC_pscores_knn_rbf_8e-1_rare.txt"
pscores.to_csv(fname, sep="\t")
print("Saved pscores")
del pscores
gc.collect()
