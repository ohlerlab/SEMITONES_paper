# python 3.6

import os
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
os.chdir("../data/processed/")  # set data directory


### LOAD CELL X PEAK ###
print("Loading peak accesibility data")
fname = "greenleaf_scATAC_peaks_filtered.mtx"
peaks = mmread(fname)  # read in the sparse matrix
peaks = peaks.tocsr()

fname = "greenleaf_scATAC_peaks_filtered.txt"
with open(fname, "r") as f:
    peaknames = [p.strip("\n") for p in f.readlines()]
f.close()


### FILTER OUT NON-PRESENT PEAKS ###
idxs = np.where(peaks.getnnz(axis=0) > 0)[0]
peaks = peaks[:, idxs]
peaks = csr_matrix(peaks)


### CREATE A HIGHLY EXPRESSED SUBSET ###
print("Creating subset of peaks found in >= 1% of cells")
idxs = np.where(peaks.getnnz(axis=0) >= (peaks.shape[0] * 0.01))[0]
sub1 = peaks[:, idxs]
sub1 = csr_matrix(sub1)

print("Saving highly expressed subset")
fname = "greenleaf_scATAC_peaks_common.mtx"
mmwrite(fname, sub1)

print("Saving common peak names")
sub1names = list(np.array(peaknames)[idxs])
fname = "greenleaf_scATAC_peaks_common_labels.txt"
with open(fname, "w") as f:
    for p in sub1names:
        f.write("{0}\n".format(p))
f.close()

### CREATE A MEDIUM EXPRESSED SUBSET ###
print("Creating subset of peaks found in 0.5%-1% of cells")
idxs = np.where((peaks.getnnz(axis=0) < (peaks.shape[0] * 0.01)) &
                (peaks.getnnz(axis=0) >= (peaks.shape[0] * 0.005)))[0]
sub2 = peaks[:, idxs]
sub2 = csr_matrix(sub2)

print("Saving lowly expressed subset")
fname = "greenleaf_scATAC_peaks_medium.mtx"
mmwrite(fname, sub2)

print("Saving rare peak names")
sub2names = list(np.array(peaknames)[idxs])
fname = "greenleaf_scATAC_peaks_medium_labels.txt"
with open(fname, "w") as f:
    for p in sub2names:
        f.write("{0}\n".format(p))
f.close()

### CREATE A LOWLY EXPRESSED SUBSET
print("Creating subset of peaks found in < 0.5% of cells")
idxs = np.where(peaks.getnnz(axis=0) < (peaks.shape[0] * 0.005))[0]
sub3 = peaks[:, idxs]
sub3 = csr_matrix(sub3)

print("Saving lowly expressed subset")
fname = "greenleaf_scATAC_peaks_rare.mtx"
mmwrite(fname, sub3)

print("Saving rare peak names")
sub3names = list(np.array(peaknames)[idxs])
fname = "greenleaf_scATAC_peaks_rare_labels.txt"
with open(fname, "w") as f:
    for p in sub3names:
        f.write("{0}\n".format(p))
f.close()
