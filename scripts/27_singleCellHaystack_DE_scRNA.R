# import libraries
library(singleCellHaystack)
library(Matrix)

# set data directory
setwd("../data/processed/")

# load the data
log_X <- readMM("greenleaf_scRNA_combined_expressed_logcnorm.mtx")
exp_genes <- read.table("greenleaf_scRNA_combined_expressed_genes.txt")
cells <- read.table("greenleaf_scRNA_combined_filtered_barcodes.txt")
log_X <- t(log_X)  # cells as rows, genes as columns
colnames(log_X) <- cells$V1
rownames(log_X) <- exp_genes$V1

MD_UMAP <- read.table("greenleaf_scRNA_combined_umap25.txt")
UMAP <- read.table("greenleaf_scRNA_combined_umap.txt")

# define whether a gene is detected
median.per.gene <- apply(log_X,1,median) # get the median read count per gene
detection <- log_X > median.per.gene # use the medians as threshold
general.detection = apply(detection, 2, sum)

# run singleCellHaystack
results <- haystack(x=MD_UMAP, detection=detection, method="highD",
                    use.advanced.sampling=general.detection)

# save results
write.table(results$results, file="greenleaf_scRNA_singleCellHaystack_DE.txt",
            sep="\t")
