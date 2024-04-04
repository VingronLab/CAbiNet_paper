library(SingleCellExperiment)
library(zellkonverter)

sce <- readH5AD("../../discussed_data/raw/brain_organoids_adata_merged_RAW.h5ad")

names(assays(sce)) <- "counts"
colnames(sce) <- sce$cell_id

# subset to triple-i protocol treated cells.
# Removes dublets and cells of unknown type. 
# remove if you want all ~100k cells.
sce <- sce[, which(sce$protocol == "Triple-i")]
sce <- sce[, which(sce$cell_type != "doublet" & sce$cell_type != "Unknown")]

# save
saveRDS(sce, file.path("../../discussed_data/raw/brain_organoids_RAW.rds"))
