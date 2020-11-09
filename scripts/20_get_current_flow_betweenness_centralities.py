# python 3.6

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx

os.chdir("../data/processed/")

files = ["greenleaf_scRNA_geneset_graphs_interaction.pkl",
         "greenleaf_scRNA_geneset_graphs_min.pkl",
         "greenleaf_scRNA_geneset_graphs_max.pkl",
         "greenleaf_scRNA_geneset_graphs_median.pkl",
         "greenleaf_scRNA_tfset_graphs_interaction.pkl",
         "greenleaf_scRNA_tfset_graphs_min.pkl",
         "greenleaf_scRNA_tfset_graphs_max.pkl",
         "greenleaf_scRNA_tfset_graphs_median.pkl"]

for file in files:

    combtyp = file.split("_")[4].split(".")[0]
    settyp = file.split("_")[2]

    with open(file, "rb") as f:
        graphs = pickle.load(f)
    f.close()

    centralities = {}
    for k, v in graphs.items():
        try:
            c = nx.current_flow_betweenness_centrality(v)
        except:
            continue
        c = {k: v for k, v in c.items() if v != 0}
        centralities[k] = c
    with open(f"greenleaf_scRNA_centralities_{settyp}_{combtyp}.pkl", "wb") as f:
        pickle.dump(centralities, f)
    f.close()
