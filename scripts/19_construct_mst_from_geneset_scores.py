# python 3.6

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from SEMITONES.enrichment_scoring import sig_interval
from SEMITONES.support_funcs import sig_dictionary

os.chdir("../data/processed/")

fname = "greenleaf_scRNA_escores_knn_rbf_8e-1.txt"
genescores = pd.read_csv(fname, sep="\t", index_col=0)

exptypes = ["interaction", "min", "max", "median"]

fname = "greenleaf_scRNA_combined_expressed_genes.txt"
with open(fname, "r") as f:
    genes = [s.strip("\n") for s in f.readlines()]
f.close()
genedict = dict(zip(range(len(genes)), genes))

for exptype in exptypes:

    fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_{exptype}_escores_all.h5"
    escores = pd.read_hdf(fname)
    fname = f"greenleaf_scRNA_gs_knn_rbf_8e-1_{exptype}_pscores_all.h5"
    pscores = pd.read_hdf(fname)
    del fname

    sigsint = sig_interval(pscores, n_sds=30)
    sigdict = sig_dictionary(escores, sigsint, sign="positive")
    del sigsint

    graphs = {}
    for c in escores.columns:
        sigs = escores.loc[sigdict[c], c]

        graphs[c] = {}
        new_index = []
        for pair in sigs.index:
            g1, g2 = genedict[pair[0]], genedict[pair[1]]
            new_index.append((g1, g2))
        sigs.index = new_index
        nodes = np.unique([i[0] for i in sigs.index] +
                          [i[1] for i in sigs.index])     
        sizes = genescores.loc[nodes, c].values
        sizednodes = [(n, {"size": s}) for n, s in zip(nodes, sizes)]
        weights = [sigs.index[i] + (sigs.iloc[i], )
                   for i in range(len(sigs))]
        g = nx.Graph()
        g.add_nodes_from(sizednodes)
        g.add_weighted_edges_from(weights)
        g = nx.algorithms.tree.mst.maximum_spanning_tree(g)
        graphs[c] = g
    fname = f"greenleaf_scRNA_geneset_graphs_{exptype}.pkl"
    with open(fname, "wb") as f:
        pickle.dump(graphs, f)
    f.close()
    del escores, pscores, graphs
