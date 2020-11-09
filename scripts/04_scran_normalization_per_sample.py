# python 3.6

import os
import gc
from scipy.io import mmread
from rpy2.robjects import numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
scran = importr("scran")

os.chdir("../data/interim/")
files = [f for f in os.listdir() if "filtered.mtx" in f]
files = [f for f in files if not "combined" in f]

frames = {}
for f in files:
    s, r = f.split("_")[2:4]
    s = "{0}_{1}".format(s, r)
    frame = mmread(f).T.todense()
    
    clusters = scran.quickCluster(frame, method=("igraph"))
    sfs = scran.calculateSumFactors(frame, clusters=clusters)
    del clusters, frame
    gc.collect()
    sfs = [sf for sf in sfs]  # make list of size factors

    # save list of size factors
    fname = "greenleaf_scRNA_{0}_scran_size_factors.txt".format(s)
    with open(fname, "w") as f:
        for sf in sfs:
            f.write("{0}\n".format(sf))
    f.close()
    print("Size factors for {0} saved".format(s))
    del sfs
    gc.collect()

numpy2ri.deactivate()
print("Done")