# python 3.6

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


### LOAD LOGC EXPRESSION ###
print("Loading logc expression")
fname = "greenleaf_scRNA_combined_expressed_logcnorm.mtx"
logcnorm = mmread(fname)
logcnorm = logcnorm.tocsr()


### LOAD GENES ###
print("Loading gene names")
fname = "greenleaf_scRNA_combined_expressed_genes.txt"
with open(fname, "r") as f:
    genes = [g.strip("\n") for g in f.readlines()]
f.close()


### REFERENCE CELLS ###
print("Obtaining reference cells")
fname = "greenleaf_scRNA_cells_selected_from_knn.txt"
with open(fname, "r") as f:
    r_cells = [int(c.strip("\n")) for c in f.readlines()]
f.close()


### SIMILARITY TO R_CELLS ###
print("Calculating similarities to the reference cells")
# load 25D UMAP
umap25 = np.load("greenleaf_scRNA_combined_umap25.npy")
metric = "rbf"  # radial basis function kernel
gamma = 8e-2  # influence radius
S = pairwise_similarities(umap25, r_cells, metric="rbf",
                          metric_params={"gamma": gamma})
del umap25
gc.collect()


### ENRICHMENT SCORING ###
print("Enrichment scoring")
scores = calculate_escores(logcnorm, query=r_cells, S=S, scale_exp=True)
scores.index = genes
fname = "greenleaf_scRNA_escores_knn_rbf_8e-2.txt"
scores.to_csv(fname, sep="\t")
print("Saved scores")
del scores
gc.collect()


### PERMUTE DATAFRAME ###
print("Permuting dataframe")
P = permute(logcnorm, n=256, axis=0, seed=42)
del logcnorm
gc.collect()


### P ENRICHMENT SCORING ###
print("Enrichment scoring over permutated dataframe")
pscores = calculate_escores(P, query=r_cells, S=S, scale_exp=True)
fname = "../interim/greenleaf_scRNA_pscores_knn_rbf_8e-2.txt"
pscores.to_csv(fname, sep="\t")
print("Saved pscores")
del pscores
gc.collect()
