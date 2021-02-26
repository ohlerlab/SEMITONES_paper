library(singleCellHaystack)

setwd("../data/simulations/")

n_cells <- c(1000, 3000, 5000, 7000, 10000)
n_clusts <- c(2, 6, 10, 14, 20)

all_times <- data.frame(it=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
for (i in seq(1, 5)) {  # for each dataset
    times <- c()
    for (j in seq(1, 10)) {  # for each iteration
        # load the expression data
        Xname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clusts[i]), "_clusts_it",
                       as.character(j), "_groups_counts.txt",
                       sep="")
        X <- read.table(Xname, sep="")

        # load 30 PCs
        PCname <- gsub("counts", "pca", Xname)
        PC30 <- read.table(PCname)

        # start timing
        Sys.sleep(60)  # sleep 1 minute
        start <- Sys.time()
        
        # get DE scores
        median.per.gene <- apply(X,1,median) # median read counts
        detection <- X > median.per.gene  # set detection cutoff
        general.detection = apply(detection, 2, sum)
        set.seed(42)
        results <- haystack(x=PC30, detection=detection, method="highD",
                            use.advanced.sampling=general.detection)
        
        # calculate times
        passed <- difftime(Sys.time(), start, unit="secs")
        passed <- as.numeric(passed)
        times <- append(times, passed)
        
        # save DE results
        savename = gsub("pca", "sCH_PCA", PCname)
        write.table(results$results, savename)
    }
    # write times to table
    colname <- paste(as.character(n_cells[i]), "_cells_",
                     as.character(n_clusts[i]), "clusts", sep="")
    all_times[,colname] <- times
}
# save timing results
write.table(all_times, "groups_sCH_PCA_times.txt", sep="\t")
