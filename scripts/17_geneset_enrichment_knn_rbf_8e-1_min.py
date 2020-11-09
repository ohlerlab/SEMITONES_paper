# python 3.6

import os
import gc
import sys
import pickle
import numpy as np
from itertools import chain
from scipy.sparse import issparse
from SEMITONES._utils import _make_generator
from SEMITONES.support_funcs import load_sparse_h5
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.enrichment_scoring import feature_set_values
from SEMITONES.enrichment_scoring import calculate_escores, permute
os.chdir("../data/processed/")
gc.enable()

i = sys.argv[1]

# load regerence cell subset
print(f"Loading reference cells ({i})")
fname = "greenleaf_scRNA_cells_selected_from_knn_subset.txt"
with open(fname, "r") as f:
    r_cells = [int(c.strip("\n")) for c in f.readlines()]
f.close()

# calculate similarities to the reference cells
print(f"Calculating similarities ({i})")
umap25 = np.load("greenleaf_scRNA_combined_umap25.npy")
S = pairwise_similarities(umap25, r_cells, metric="rbf",
                          metric_params={"gamma": 8e-1})
del umap25
gc.collect()

# load the expression data
print(f"Loading expression data ({i})")
fname = "greenleaf_scRNA_combined_expressed_logcnorm.h5"
data = load_sparse_h5("scRNA", fname)
if issparse(data) and (data.getformat() not in ["csr", "csc"]):
    data = data.tocsr()

# make set dictionary working directory
os.chdir("set_scores/")
    
# load the gene sets
print(f"Loading gene sets ({i})")
fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_idxs_{i}.pkl"
with open(fname, "rb") as f:
    genesets = pickle.load(f)
f.close()
genesets = _make_generator(genesets)

# create gene set expression frame
print(f"Construction gene set expression frame ({i})")
gsexp = feature_set_values(data, genesets, combtype="min")
del data
gc.collect()

# calculate enrichment scores
print(f"Calculating escores ({i})")
escores = calculate_escores(gsexp, query=r_cells, S=S, scale_exp=True)
fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_min_escores_{i}.txt"
np.savetxt(fname, escores.values, delimiter="\t")
del escores
gc.collect()

# permute gsexp
print(f"Permuting gene set expression frame ({i})")
P = permute(gsexp, n=256, axis=0, seed=42)
del gsexp
gc.collect()

# calculate pscores
print(f"Calculating pscores ({i})")
pscores = calculate_escores(P, query=r_cells, S=S, scale_exp=True)
fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_min_pscores_{i}.txt"
np.savetxt(fname, pscores.values, delimiter="\t")
