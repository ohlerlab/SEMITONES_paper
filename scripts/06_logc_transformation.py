# python 3.6

import os
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
os.chdir("../data/processed/")  # set data directory


### LOAD CELL X GENE EXPRESSION ###
print("Loading expression data")
fname = "greenleaf_scRNA_combined_norm.mtx"
norm = mmread(fname)  # read in the sparse matrix
norm = norm.tocsr()


### ONLY KEEP EXPRESSED GENES ###
print("Removing non-expressed genes")
idxs = np.where(norm.getnnz(axis=0) > 0)[0]
norm = norm[:, idxs]

print("Saving gene subset list")
fname = "../interim/greenleaf_scRNA_combined_genes.txt"
with open(fname, "r") as f:
    genes = [g.strip("\n") for g in f.readlines()]
f.close()

genes = list(np.array(genes)[idxs])

fname = "greenleaf_scRNA_combined_expressed_genes.txt"
with open(fname, "w") as f:
    for g in genes:
        f.write("{0}\n".format(g))
f.close()


### LOG-TRANSFORMATION ###
print("Log-transforming expression data Lun-way")

fname = "../interim/greenleaf_scRNA_combined_scran_size_factors.txt"
with open(fname, "r") as f:
    sfs = [float(s.strip("\n")) for s in f.readlines()]  # get size factors
f.close()

smin, smax = np.percentile(sfs, q=[5, 95])
delta = np.abs([smin ** -1 - smax ** -1])
pseudo = np.max([1, delta])  # rho is 1 so not specified

logcnorm = np.log(norm + pseudo)
logcnorm = csr_matrix(logcnorm)  # make sure it is the right format

print("Saving log transformed and filtered expression data")
mmwrite("greenleaf_scRNA_combined_expressed_logcnorm.mtx", logcnorm)
