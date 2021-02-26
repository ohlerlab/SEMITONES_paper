library(Seurat)
setwd("../data/simulations/")

n_cells = c(1000, 3000, 5000, 7000, 10000)
n_clust = c(2, 6, 10, 14, 20)

# grouped data
for (i in seq(1, 5)) {
    for (j in seq(1, 10)) {
        fname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clust[i]), "_clusts_it",
                       as.character(j), "_groups_counts.txt", sep="")
        
        print(paste("Processing ", fname))
        
        counts <- read.table(fname)
        
        obj <- CreateSeuratObject(counts, project="sim", verbose=FALSE,
                                  min.cells=1)  # gene must be expressed
        obj <- NormalizeData(obj, verbose=FALSE)
        obj <- ScaleData(obj, verbose=FALSE)
       
        # save the normalized and scaled data
        fname <- gsub("counts", "scaledata", fname)
        write.table(obj@assays$RNA@scale.data, fname)                    
                    
        # save the top 30 PCs
        obj <- RunPCA(obj, features=rownames(counts), npcs=30, verbose=FALSE)
        fname <- gsub("scaledata", "pca", fname)
        write.table(Embeddings(obj, reduction="pca"), fname)

        obj <- FindNeighbors(obj, verbose=FALSE)
        obj <- FindClusters(obj, verbose=FALSE)

        # save the 2D UMAP
        obj <- RunUMAP(obj, dims=1:30, verbose=FALSE)
        fname <- gsub("pca", "umap", fname)
        write.table(Embeddings(obj, reduction="umap"), fname)
        
        # save the 10D UMAP
        obj <- RunUMAP(obj, dims=1:30, max.dim=10, verbose=FALSE)
        fname <- gsub("umap", "umap10", fname)
        write.table(Embeddings(obj, reduction="umap"), fname)

        # save the Seurat object
        fname <- gsub("umap10", "seuratobj", fname)
        fname <- gsub(".txt", ".rds", fname)
        saveRDS(obj, fname)
    }
}

# path data
for (i in seq(1, 5)) {
    for (j in seq(1, 10)) {
        fname <- paste(as.character(n_cells[i]), "_cells_",
                       as.character(n_clust[i]), "_clusts_it",
                       as.character(j), "_paths_counts.txt", sep="")
        
        print(paste("Processing ", fname))
        
        counts <- read.table(fname)
        
        obj <- CreateSeuratObject(counts, project="sim", verbose=FALSE,
                                  min.cells=1)  # gene must be expressed
        obj <- NormalizeData(obj, verbose=FALSE)
        obj <- ScaleData(obj, verbose=FALSE)

        # save the normalized and scaled data
        fname <- gsub("counts", "scaledata", fname)
        write.table(obj@assays$RNA@scale.data, fname)                    

        # save the top 30 PCs
        obj <- RunPCA(obj, features=rownames(counts), npcs=30, verbose=FALSE)
        fname <- gsub("scaledata", "pca", fname)
        write.table(Embeddings(obj, reduction="pca"), fname)

        obj <- FindNeighbors(obj, verbose=FALSE)
        obj <- FindClusters(obj, verbose=FALSE)

        # save the 2D UMAP
        obj <- RunUMAP(obj, dims=1:30, verbose=FALSE)
        fname <- gsub("pca", "umap", fname)
        write.table(Embeddings(obj, reduction="umap"), fname)

        # save the 10D UMAP
        obj <- RunUMAP(obj, dims=1:30, max.dim=10, verbose=FALSE)
        fname <- gsub("umap", "umap10", fname)
        write.table(Embeddings(obj, reduction="umap"), fname)

        # save the Seurat object
        fname <- gsub("umap10", "seuratobj", fname)
        fname <- gsub(".txt", ".rds", fname)
        saveRDS(obj, fname)
    }
}