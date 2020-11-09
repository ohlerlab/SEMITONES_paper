# python 3.6

import os
import sys
import timeit as t
import numpy as np
import pandas as pd
from functools import partial
from scipy.sparse import random, vstack
from SEMITONES.enrichment_scoring import calculate_escores
os.chdir("../data/processed/")

# write timing function
def timeit(X, S, r_cells, ncpu, optim_over):
    
    # define the function to time test
    f = partial(calculate_escores,  # function to compute
                X, query=r_cells, S=S,  # inpute values
                ncpu=ncpu, scale_exp=True, optim_over=optim_over,
                make_copy=False)  # specs

    time = t.repeat(f, repeat=3, number=1)  # time a single run thrice

    return time

# define parameters using passed parameters
s = float(sys.argv[1])  # sparsity character
n_cells = int(sys.argv[2])  # list of number of cells in total
n_rcells = sys.argv[3]  # number of reference cells
n_rcells = map(int, n_rcells.strip("[]").split(","))
n_rcells = np.array(list(n_rcells))
n_features = sys.argv[4]  # list of number of features
n_features = map(int, n_features.strip("[]").split(","))
n_features = np.array(list(n_features))
max_ncpu = int(sys.argv[5])  # maximum number of CPUs to use
optim_over = sys.argv[6]
outfile = sys.argv[7]

# define dataframes to subset
nfmax = int(n_features.max())
if (n_cells > 400000):
    if n_cells % 2 != 0:
        print("Please provide an even number o f cells (n_cells)")
    n = int(n_cells / 20)
    X = random(n, nfmax, s, "csr")
    X = vstack([X, X, X, X, X, X, X, X, X, X,
                X, X, X, X, X, X, X, X, X, X])
elif (nfmax < 250_000) & (s < 1):
    X = random(n_cells, nfmax, s, "csr")
elif (nfmax < 400_000) & (s < 1):
    # check if the number of cells is even
    if n_cells % 2 != 0:
        print("Please provide an even number of cells (n_cells)")
    n = int(n_cells / 2)
    X = random(n, nfmax, s, "csr")
    X = vstack([X, X])

else:
    # check if the number of cells if even
    if n_cells % 2 != 0:
        print("Please provide an even number of cells (n_cells)")
    n = int(n_cells / 4)
    X = random(n, nfmax, s, "csr")
    X = vstack([X, X, X, X])
X.eliminate_zeros()
print(X.shape)

# do this function for a certain number of reference cells and features
times = {}
for n_r in n_rcells:
    n_r = int(n_r)
    for n_f in n_features:
        n_f = int(n_f)

        # check is ncpu <= n_features if optim_over == "rows"
        if (optim_over == "rows") and (max_ncpu > n_r):
            ncpu = n_r
        elif (optim_over == "cols") and (max_ncpu > n_f):
            ncpu = n_f
        else:
            ncpu = max_ncpu

        # construct a similarity matrix
        S = np.random.uniform(0.0, 1.0, (n_cells, n_r))

        print(f"Timing for {n_r} reference cells and {n_f} features." + 
              f"and {ncpu} CPU cores.")

        r_cells = list(range(n_r))
        times[f"{n_r}_{n_f}"] = timeit(X[:, 0:n_f], S, r_cells,
                                       ncpu, optim_over)
times = pd.DataFrame(times)
times.to_csv(outfile, sep="\t")
