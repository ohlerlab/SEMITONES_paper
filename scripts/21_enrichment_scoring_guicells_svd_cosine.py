# python 3.6

import os
import gc
import numpy as np
import pandas as pd
from scipy.io import mmread
from SEMITONES.support_funcs import pairwise_similarities, load_sparse_h5
from SEMITONES.enrichment_scoring import calculate_escores, permute

os.chdir("../data/processed/")

fname = "greenleaf_scRNA_combined_expressed_logcnorm.h5"
X = load_sparse_h5("scRNA", fname)

SVD = mmread("greenleaf_scRNA_combined_svd50.mtx")

fname = "greenleaf_scRNA_cells_selected_from_gui.txt"
with open(fname, "r") as f:
    rcells = [int(c.strip("\n")) for c in f.readlines()]
f.close()

fname = "greenleaf_scRNA_combined_expressed_genes.txt"
with open(fname, "r") as f:
    genes = [g.strip("\n") for g in f.readlines()]
f.close()

S = pairwise_similarities(SVD, rcells, metric="cosine")

escores = calculate_escores(X, rcells, S=S, scale_exp=True)
escores.index = genes
escores.to_csv("greenleaf_scRNA_escores_knn_cosineoversvd.txt", sep="\t")
del escores

P = permute(X, n=256, seed=42)
del X

pscores = calculate_escores(P, rcells, S=S, scale_exp=True)
pscores.to_csv("../interim/greenleaf_scRNA_pscores_knn_cosineoversvd.txt", sep="\t")
