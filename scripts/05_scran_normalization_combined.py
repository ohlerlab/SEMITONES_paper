# python 3.6

import os
import gc
from scipy.io import mmread
from rpy2.robjects import numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
scran = importr("scran")

os.chdir("/fast/AG_Ohler/anna/hcv1/data/interim/")
fname = "greenleaf_scRNA_combined_filtered.mtx"
frame = mmread(fname).todense().T

clusters = scran.quickCluster(frame, method=("igraph"))
sfs = scran.calculateSumFactors(frame, clusters=clusters)
sfs = [sf for sf in sfs]  # get list of size factors

# save list of size factors
fname = "greenleaf_scRNA_combined_scran_size_factors.txt"
with open(fname, "w") as f:
    for sf in sfs:
        f.write("{0}\n".format(sf))
f.close()
numpy2ri.deactivate()
print("Done")