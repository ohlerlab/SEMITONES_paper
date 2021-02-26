library(Seurat)

setwd("../data/simulations/")

n_cells <- c(1000, 3000, 5000, 7000, 10000)
n_clusts <- c(2, 6, 10, 14, 20)

all_times <- data.frame(it=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
for (i in seq(1, 5)) {  # for each dataset
    times <- c()
    for (j in seq(1, 10)) {  # for each iteration
        # read the seurat object
        Xname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clusts[i]), "_clusts_it",
                       as.character(j), "_paths_seuratobj.rds",
                       sep="")
        X <- readRDS(Xname)

        # use the clusters as identity
        Idents(X) <- X@meta.data$RNA_snn_res.0.8

        # start the timer
        Sys.sleep(60)  # sleep 1 minute
        start <- Sys.time()

        # get the DE scores
        DE <- FindAllMarkers(X, return.thresh=1)

        # calculate the passed times
        passed <- difftime(Sys.time(), start, unit="secs")
        passed <- as.numeric(passed)
        times <- append(times, passed)

        # save the DE scores
        savename = gsub("seuratobj.rds", "Wilcoxon_filtered.txt", Xname)
        write.table(DE, savename)
    }
    # write the timing results to the table
    colname <- paste(as.character(n_cells[i]), "_cells_",
                     as.character(n_clusts[i]), "clusts", sep="")
    all_times[,colname] <- times
}
# save the timing results
write.table(all_times, "paths_Wilcoxon_filtered_times.txt", sep="\t")