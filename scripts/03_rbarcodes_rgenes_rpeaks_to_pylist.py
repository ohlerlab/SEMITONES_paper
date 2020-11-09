# python 3.6

import os
import pandas as pd

os.chdir("../data/interim")

# get all barcode file names
bfnames = [f for f in os.listdir() if "barcodes" in f]
for f in bfnames:
    barcodes = pd.read_csv(f).V1  # get the barcodes as series
    barcodes = [f.split(" ")[1].strip('"') for f in barcodes]  # make list
    # save barcodes as list in simple .txt file
    with open(f, "w") as f:
        for b in barcodes:
            f.write("{0}\n".format(b))
    f.close()

# get all gene name files
gfnames = [f for f in os.listdir() if "genes" in f]
for f in gfnames:
    genes = pd.read_csv(f).iloc[:, 0]  # get gene names
    genes = [f.split(" ")[1].strip('"') for f in genes]  # make list
    # save genes as list in simple .txt files
    with open(f, "w") as f:
        for g in genes:
            f.write("{0}\n".format(g))
    f.close()

# get all peak name files
pfnames = [f for f in os.listdir() if "peaks" in f]
for f in pfnames:
    peaks = pd.read_csv(f).iloc[:, 0]  # get peak names
    peaks = [f.split(" ")[1].strip('"') for f in peaks]  # make list
    # save list of peak names as simple .txt files
    with open(f, "w") as f:
        for p in peaks:
            f.write("{0}\n".format(p))
    f.close()
