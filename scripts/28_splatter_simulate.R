library(splatter)
setwd("../data/simulations/")

n_cells = c(1000, 3000, 5000, 7000, 10000)
n_clust = c(2, 6, 10, 14, 20)
seeds = c(32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384)

# create cluster/"groups" simulations
for (i in seq(1, 5)) {  # idx number cell and cluster
    for (j in seq(1, 10)) {  # idx seed
        params <- newSplatParams(seed=seeds[j])
        group.prob <- rep(1/n_clust[i], n_clust[i])  # prop in cluster
        dataset <- splatSimulate(params, group.prob=group.prob,
                                 method="groups", batchCells=n_cells[i],
                                 seed=seeds[j])
        # save the counts
        fname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clust[i]), "_clusts_it",
                       as.character(j), "_groups_counts.txt", sep="")
        write.table(assays(dataset)$counts, fname)
        # save the column metadata
        fname <- gsub("counts.txt", "coldata.txt", fname)
        write.table(colData(dataset), fname)
        # save the row metadata
        fname <- gsub("coldata.txt", "rowdata.txt", fname)
        write.table(rowData(dataset), fname)
    }
}

# create trajectory/"path" simulations
for (i in seq(1, 5)) {  #
    for (j in seq(1, 10)) {
        params <- newSplatParams(seed=seeds[j])
        group.prob <- rep(1/n_clust[i], n_clust[i])  # prop in trajectory
        dataset <- splatSimulate(params, group.prob=group.prob,
                                 method="paths", batchCells=n_cells[i],
                                 seed=seeds[j])
        # save the counts
        fname <- paste(as.character(n_cells[i]), "_cells_",
                      as.character(n_clust[i]), "_clusts_it",
                      as.character(j), "_paths_counts.txt", sep="")
        write.table(assays(dataset)$counts, fname)
        # save the column metadata
        fname <- gsub("counts.txt", "coldata.txt", fname)
        write.table(colData(dataset), fname)
        # save the row metadata
        fname <- gsub("coldata.txt", "rowdata.txt", fname)
        write.table(rowData(dataset), fname)
    }
}
