# to use this script make sure the raw data is in /data/raw
# R 3.6.2

library(Matrix)
library(SummarizedExperiment)

setwd("../data/")

fnames <- list.files(path="raw/", pattern="*scRNA*")

save_mtx_barcodes_genes <- function(filename) {
    splitname <- strsplit(filename, "_")[[1]]
    origin <- splitname[3]
    sample <- strsplit(splitname[4], "[.]")[[1]][1]

    # load the RDS object with the dataframe
    readname <- paste0("raw/", filename)
    df <- readRDS(readname)

    # save the transposed dataframe
    writename <- paste0("interim/greenleaf_scRNA_", origin, "_",
                        sample, ".mtx")
    writeMM(t(df), writename)

    # save the barcodes
    bfname <- paste0("interim/greenleaf_scRNA_", origin, "_",
                     sample, "_barcodes.txt")
    write.table(as.matrix(colnames(df)), file=bfname)

    # save the gene names
    gfname <- paste0("interim/greenleaf_scRNA_", origin, "_",
                     sample, "_genes.txt")
    write.table(as.matrix(rownames(df)), file=gfname)
}

lapply(fnames, save_mtx_barcodes_genes)