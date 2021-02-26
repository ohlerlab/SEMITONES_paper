# import libraries
library(Seurat)
library(Matrix)

# set data location
setwd("../data/processed/")

# load normalized count data
print("Read count data")
X <- readMM("greenleaf_scRNA_combined_norm.mtx")
genes <- read.table("greenleaf_scRNA_combined_genes.txt")  # genes
cells <- read.table("greenleaf_scRNA_combined_filtered_barcodes.txt")  # cells
X <- t(X)  # rows are genes, columns are cells
rownames(X) <- genes$V1
colnames(X) <- cells$V1

# load log-transformed count data (because we used an alternative pseudocount)
print("Read log-transformed data")
log_X <- readMM("greenleaf_scRNA_combined_expressed_logcnorm.mtx")
exp_genes <- read.table("greenleaf_scRNA_combined_expressed_genes.txt")
log_X <- t(log_X)  # rows are genes, columns are cells
rownames(log_X) <- exp_genes$V1
colnames(log_X) <- cells$V1

# filter non-log-transformed count data for expressed genes
print("Filter genes")
X <- X[exp_genes$V1,]

# load metadata
print("Load metadata")
metadata <- read.table("../external/greenleaf_scRNA_cell_metadata.txt")
newcols <- paste(metadata$Group, metadata$Barcode, sep=":")
rownames(metadata) <- newcols
metadata <- metadata[colnames(X),]

# create Seurat object
print("Create Seurat object")
obj <- CreateSeuratObject(counts=X, project="DE_testing", assay="RNA",
                          meta.data=metadata, min.cells=0, min.features=0)

# add log-transformed data into data slot
print("Add log-transformed data")
log_X <- as(object=log_X, Class='dgCMatrix')
obj <- SetAssayData(object=obj, assay="RNA", new.data=log_X,
                    slot="data")

# set clusters as cell identity
print("Set cell identity")
Idents(obj) <- obj@meta.data$Clusters

# run FindAllMarkers for all genes and all clusters
print("Identify markers")
DE <- FindAllMarkers(obj, assay="RNA", slot="data", return.thresh=1)

# save DE genes to text df
write.table(DE, file="greenleaf_scRNA_Seurat_DE_filtered.txt", sep="\t")
