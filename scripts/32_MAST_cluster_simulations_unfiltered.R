library(Seurat)

setwd("../data/simulations/")

n_cells <- c(1000, 3000, 5000, 7000, 10000)
n_clusts <- c(2, 6, 10, 14, 20)

all_times <- data.frame(it=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
for (i in seq(1, 5)) {  # for each dataset
    times <- c()
    for (j in seq(1, 10)) {  # for each iteration
        # read Seurat object
        Xname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clusts[i]), "_clusts_it",
                       as.character(j), "_groups_seuratobj.rds",
                       sep="")
        X <- readRDS(Xname)

        # use clusters as identity
        Idents(X) <- X@meta.data$RNA_snn_res.0.8

        # find marker genes
        DE <- FindAllMarkers(X, test.use="MAST", return.thresh=1,
                             min.pct=0, logfc.threshold=0)
        
        # save teh results
        savename = gsub("seuratobj.rds", "MAST_filtered.txt", Xname)
        write.table(DE, savename)
    }
}
