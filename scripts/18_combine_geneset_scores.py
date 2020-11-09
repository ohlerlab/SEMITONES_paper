# python 3.6

import os
import gc
import pickle
import numpy as np
import pandas as pd
from SEMITONES.support_funcs import save_sparse_h5

os.chdir("../data/processed/set_scores")

def combine_all_scores(settype, combtype, scoretype, n):
    """ Combines dataframes and index files into one
    settype: str
        The set type to indicate in the file
    combtype: str
        The combination type as set in feature_set_values()
    scoretype: str
        Whether the scores a real scores or the null distribution:
        "pscores" for the null, "escores" for actual scores
    n: int
        How many files need to be combined (as indexed in their name)
       """
    idxs = []
    scores = []
    for i in range(1, n + 1):
        fname = f"greenleaf_scRNA_{settype}_knn_rbf_8e-1_idxs_{i}.pkl"
        with open(fname, "rb") as f:
            idxs.extend(pickle.load(f))
        f.close()

        fname = f"greenleaf_scRNA_{settype}_knn_rbf_8e-1_{combtype}_{scoretype}_{i}.txt"
        scores.append(np.loadtxt(fname))

    scores = np.vstack(scores)
    scores = pd.DataFrame(scores, index=idxs)
    del idxs
    gc.collect()
    fname = f"../greenleaf_scRNA_cells_selected_from_knn_subset.txt"
    with open(fname, "r") as f:
        columns = [s.strip("\n") for s in f.readlines()]
    f.close()
    scores.columns = columns
    fname = f"../greenleaf_scRNA_{settype}_knn_rbf_8e-1_{combtype}_{scoretype}_all.h5"
    scores.to_hdf(fname, key="scores", mode="w")


for settype, n in zip(["gs"], [67]):
    for combtype in ["interaction", "min", "max", "median"]:
        for scoretype in ["escores", "pscores"]:
            combine_all_scores(settype, combtype, scoretype, n)
