library(Seurat)
library(SeuratWrappers)
library(monocle3)

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

        # import into monocle3 and process
        X <- as.cell_data_set(X)  # from SeuratWrappers
        X <- cluster_cells(X)
        
        # start time
        Sys.sleep(60)
        start <- Sys.time()
        
        # compute graph scores
        X <- learn_graph(X)
        graphscores <- graph_test(X)
        
        # calculate time
        passed <- difftime(Sys.time(), start, unit="secs")
        passed <- as.numeric(passed)
        times <- append(times, passed)

        # save results
        graphname <- gsub("seuratobj.rds", "monocle3_graphs.txt", Xname)
        write.table(graphscores, graphname)
    }
    # save times to table
    colname <- paste(as.character(n_cells[i]), "_cells_",
                     as.character(n_clusts[i]), "clusts", sep="")
    all_times[,colname] <- times
}

# save timings
write.table(all_times, "groups_monocle3_graphs_times.txt", sep="\t")
