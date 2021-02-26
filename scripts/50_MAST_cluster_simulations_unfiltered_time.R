library(Seurat)

setwd("../data/simulations/")

n_cells <- c(1000, 3000, 5000, 7000, 10000)
n_clusts <- c(2, 6, 10, 14, 20)

# because this takes so long we loop over iterations first and save each it
all_times <- data.frame(it=c(1, 2, 3, 4, 5))
for (j in seq(1, 10)) {  # for each iteration
    times <- c()
    for (i in seq(1, 5)) {  # for each dataset
        # read the Seurat object
        Xname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clusts[i]), "_clusts_it",
                       as.character(j), "_groups_seuratobj.rds",
                       sep="")
        X <- readRDS(Xname)

        # set clusters as the identity
        Idents(X) <- X@meta.data$RNA_snn_res.0.8

        # start the timer
        Sys.sleep(60)  # sleep 1 minute
        start <- Sys.time()

        #  find DE genes
        DE <- FindAllMarkers(X, test.use="MAST", return.thresh=1,
                             min.pct=0, logfc.threshold=0)

        # get the passed time
        passed <- difftime(Sys.time(), start, unit="secs")
        passed <- as.numeric(passed)
        times <- append(times, passed)
    }
    # save the times to a dataframe
    colname <- paste("it", as.character(j), sep="")
    all_times[,colname] <- times
    # save this iteration to file
    filename <- paste("groups_MAST_unfiltered_times_it",
                      as.character(j), ".txt", sep="")
    write.table(all_times, filename, sep="\t")
}
