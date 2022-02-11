import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc

from adata import AnnData
from scipy.io import mmread
from scipy.sparse import csr_matrix
from sklearn.neighbors import KDTree
from SEMITONES.enrichment_scoring import pvals_per_gene
from SEMITONES.tfidf import TFIDF
from sklearn.decomposition import TruncatedSVD
from umap import UMAP


def density_index(X, k=None):

    """Taken from https://www.biorxiv.org/content/10.1101/2020.10.07.330563v1
    Inspired by the root-mean-squared distance drms between pairs of cells
    i and j can be calculated as."""
    
    k = 5 if k is None else k

    km = np.mean(KDTree(X).query(X, k=k)[0])
    fronorm = np.linalg.norm(X, ord="fro")
    N = X.shape[0]

    return np.sqrt(2 / N) * (fronorm / km)


# get genes
os.chdir("../../data/processed/")
with open("greenleaf_scRNA_combined_expressed_genes.txt", "r") as f:
    genes = [g.strip("\n") for g in f.readlines()]
f.close()

# load the normalized raw data
fname = "greenleaf_scRNA_combined_norm.mtx"
X = mmread("greenleaf_scRNA_combined_norm.mtx")
X = csr_matrix(X)  # make csr
# get gene names
with open("greenleaf_scRNA_combined_genes.txt", "r") as f:
    genes = [g.strip("\n") for g in f.readlines()]
f.close()
# get expressed genes
with open("greenleaf_scRNA_combined_expressed_genes.txt", "r") as f:
    exp_genes = [g.strip("\n") for g in f.readlines()]
f.close()
exp_idx = [genes.index(i) for i in exp_genes]  # get indices of expressed genes
X = X[:, exp_idx]  # subset expressed genes
genes = exp_genes  # genes is expressed genes

# set directory for results
os.chdir("../../results/evaluate_feature_selection")

# use scanpy to get seurat v3 flavoured features
adata = AnnData(X)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=15000)
seurat = adata.var

# get the KS-test genes
semitones = pvals_per_gene(escores, pscores).sort_values(by="q")

# use scanpy to get seurat MVP features
adata = AnnData(X)
sc.pp.log1p(adata)

n_genes = [100, 300, 500, 1000, 3000, 5000, 10000]

DI = {}
for n in n_genes:
    DI[n] = {}
    
    # get SEMITONES KS-test indices
    sem_genes = semitones.sort_values(by="q").index[:n]
    sem_idx = [genes.index(i) for i in sem_genes]

    # get Seurat v3 indices
    seur3_idx = seurat.sort_values(by="highly_variable_rank").index[:n]
    
    # get Seurat indices
    seur_idx = sc.pp.highly_variable_genes(adata, flavor="seurat",
                                           n_top_genes=n)
    seur_idx = np.where(seur_idx.highly_variable.values)[0]
    
    # get SEMITONES ranked indices (from notebook)
    with open(f"semitones_ranked_feature_selection_{n}.pkl", "rb") as f:
        sem_rank = pickle.load(f)
    f.close()
    sem_rank_idx = [genes.index(i) for i in sem_rank]

    # name methods and make list of indices
    methods = ["semitones", "seurat", "seuratv3", "sem_rank"]
    indices = [sem_idx, seur_idx, seur3_idx, sem_rank_idx]

    # for each method name and the corresponding indices
    for m, gset in zip(methods, indices):
        # do LSI to get 50 components
        tfidf = TFIDF(X[:, gset])
        svd = TruncatedSVD(n_components=50).fit_transform(tfidf)
        np.save(f"{m}_svd50_{n}genes.npy", svd)  # save the SVD
        
        # get a UMAP with same parameters
        umap = UMAP(n_components=2, n_neighbors=30,
                    min_dist=0.3, random_state=42).fit_transform(svd)
        np.save(f"{m}_umap_{n}genes.npy", umap)  # save UMAP

        # compute density indices
        DI[n][m] = density_index(svd)

# save the density indices
with open("density_indices_RNA.pkl", "wb") as f:
    pickle.dump(DI, f)
f.close()
