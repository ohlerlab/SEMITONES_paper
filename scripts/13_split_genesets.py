# python 3.6

import os
import pickle
import numpy as np
from itertools import chain
from SEMITONES._utils import _make_generator, _chunk_generator
os.chdir("../data/processed/")

# import gene sets
fname = "greenleaf_scRNA_knn_subset_sd25_8e-1_genesets.pkl"
with open(fname, "rb") as f:
    genesets = pickle.load(f)
f.close()
genesets = chain.from_iterable([*genesets.values()])
genesets = np.unique([str(s) for s in genesets])
genesets = [(int(s.split(",")[0].strip("(")),
             int(s.split(",")[1].strip(")"))) for s in genesets]
genesets = _make_generator(genesets)

# split into subsets
n = 5_000  # number of sets per sample
subsets = _chunk_generator(genesets, size=n)
del genesets

# save to idx files
os.chdir("set_scores/")
for i, subset in enumerate(subsets, 1):
    fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_idxs_{i}.pkl"
    with open(fname, "wb") as f:
        pickle.dump(list(subset), f)
    f.close()
