### submitted to an HPC compute node requesting a single 3.20GHz core###

import os
import time
import numpy as np
import pandas as pd

from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.enrichment_scoring import calculate_escores
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.enrichment_scoring import permute


from sklearn.metrics.pairwise import pairwise_kernels

os.chdir("../data/simulations/")

n_cells = [1000, 3000, 5000, 7000, 10000]
n_clusts = [2, 6, 10, 14, 20]

times = {}
times_with_p = {}
for i in range(len(n_cells)):  # for each dataset
    times[i] = {}
    times_with_p[i] = {}
    for j in range(10):  # for each iteration
        prefix = f"{n_cells[i]}_cells_{n_clusts[i]}_clusts_it{j+1}"

        # read the normalized expression data
        Xname = prefix + "_paths_scaledata.txt"
        X = pd.read_csv(Xname, sep=" ")
        X = X.replace(",", ".", regex=True)
        X = X.astype(float)

        # load the UMAP
        UMAPname = prefix + "_paths_umap10.txt"
        UMAP10 = pd.read_csv(UMAPname, sep=" ")
        UMAP10 = UMAP10.replace(",", ".", regex=True)
        UMAP10 = UMAP10.astype(float)
        UMAP10 = UMAP10.values

        # prepare the expression data for SEMITONES
        X = X.T  # genes as columns/features
        genes, cells = X.columns, X.index  # save names
        X = X.values  # make numpy array

        # get the starting cell
        n_ret = int(X.shape[1] * 5e-3)  # return .5% of cells
        start = np.argmin(pairwise_kernels(UMAP10, metric="rbf",
                                           gamma=.1).mean(1))

        # start the timer
        time_start = time.time()

        # get reference cells
        r_cells = from_knn_dist(UMAP10, start=start,
                                n_ret=n_ret, metric="rbf",
                                metric_params={"gamma": .1})

        # calculate similarity to the reference cells
        S = pairwise_similarities(UMAP10, query=r_cells, metric="rbf",
                                  metric_params={"gamma": .1})

        # calculate enrichment scores using 1 CPU
        escores = calculate_escores(X, query=r_cells, S=S, ncpu=1,
                                    scale_exp=False)
        escores.columns = cells[r_cells]  # columns are reference cells
        escores.index = genes  # rows are genes

        # get the passed time
        times[i][j] = (time.time() - time_start)

        # compute the permutation null distribution
        P = permute(X)
        pscores = calculate_escores(P, query=r_cells, S=S, ncpu=1,
                                    scale_exp=False)
        pscores.columns = cells[r_cells]
        pscores.index = genes

        # get times with permutation scoring
        times_with_p[i][j] = (time.time() - time_start)
        
# save the times for enrichment scoring
pd.DataFrame(times).to_csv("paths_SEMITONES_UMAP_rbf_times.txt", sep="\t")
# save the permutation enrichment scores
pd.DataFrame(times_with_p).to_csv("paths_SEMITONES_UMAP_rbf_timesp.txt",
                                  sep="\t")
