# to use this script make sure the raw data is in /data/raw

library(Matrix)
library(SummarizedExperiment)

setwd("../data/raw/")

fname = "greenleaf_scATAC_cell_by_peak_processed.rds"

# load the RDS SummarizedExperiment
df <- readRDS(fname)
df <- assay(df, "counts")  # get the count dataframe

# save the count dataframe
writename <- "../interim/greenleaf_scATAC_combined.mtx"
writeMM(t(df), writename)

# save the barcodes
bfname <- "../interim/greenleaf_scATAC_combined_barcodes.txt"
write.table(as.matrix(colnames(df)), file=bfname)

# save the peak names
gfname <- "../interim/greenleaf_scATAC_combined_peaks.txt"
write.table(as.matrix(rownames(df)), file=gfname)
