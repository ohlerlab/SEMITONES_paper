# python 3.6

import os
from scipy.io import mmread
from SEMITONES.support_funcs import save_sparse_h5

os.chdir("../data/processed/")
scRNA = mmread("greenleaf_scRNA_combined_expressed_logcnorm.mtx")
save_sparse_h5(scRNA, "scRNA", "greenleaf_scRNA_combined_expressed_lognorm.h5")

os.chdir("../data/interim/")
scATAC = mmread("greenleaf_scATAC_peaks_filtered.mtx")
scATAC = scATAC.tocsr()
save_sparse_h5(scATAC, "scATAC", "greenleaf_scATAC_peaks_filtered.h5")
