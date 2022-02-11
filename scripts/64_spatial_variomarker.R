library(Seurat)  # load Seurat

# get the data rds
data = readRDS("../data/interim/10X_spatial_musbrain_ant1_norm.rds")

# run spatial selection method for all HVGs
data <- FindSpatiallyVariableFeatures(
    data, assay="SCT",
    features = VariableFeatures(data)[1:2449],
    selection.method="markvariogram"
)

# get signficant features
var <- SpatiallyVariableFeatures(data)

# write to list
lapply(var, write, "../results/spatially_variable_features_Seurat.txt",
       append=TRUE)
