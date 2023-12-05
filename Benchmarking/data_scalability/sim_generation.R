library(splatter)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(scran)
library(scRNAseq)
library(TENxPBMCData)
library(Seurat)

set.seed(12345)


# # Uncomment for splatter data based on Zeisel Brain data.
# dataset = "zeisel"

# Uncomment for splatter data based on Pbcm3k data.
dataset = "pbmc3k"

# maindir <- paste0("./", dataset, "/")
maindir <- '../../Data/sim_data/raw/data_scalability'
dir.create(maindir)

if(dataset == "zeisel"){

  sce <- ZeiselBrainData()
  clust.sce<- quickCluster(sce) 
  sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)
  sce <- logNormCounts(sce)
  sce <- runUMAP(sce)
  # plotUMAP(sce, colour_by = "level1class")

} else if (dataset == "pbmc3k"){

  sce <- TENxPBMCData(dataset = "pbmc3k")
  rownames(sce) <- rowData(sce)$Symbol_TENx
  colnames(sce) <- colData(sce)$Barcode

  # For pre-processing we used Seurat.
  pbmc <- CreateSeuratObject(counts = as.matrix(counts(sce)),
                            assay = "RNA",
                            project = "pbmc3k",
                            min.cells = 3,
                            min.features = 200,
                            meta.data = as.data.frame(colData(sce)))


  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

  # Filter data
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  no_zeros_rows <- rowSums(pbmc, slot = "counts") > 0
  pbmc <- pbmc[no_zeros_rows,]

  # Normalization
  pbmc <- NormalizeData(pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        verbose = FALSE)

  pbmc <- FindVariableFeatures(pbmc,
                               nfeatures = 2000,
                               verbose = FALSE)

  # Scaling
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc,
                    features = all.genes,
                    verbose = FALSE)

  # Run PCA
  pbmc <- RunPCA(pbmc,
                 features = VariableFeatures(object = pbmc),
                 verbose = FALSE)

  # Cell clustering
  pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
  pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

  pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
  # DimPlot(pbmc, reduction = "umap", label = FALSE, pt.size = 0.5)

  new.cluster.ids <- c("naive_CD4_Tcell",
                       "CD14_monocyte", 
                       "memory_CD4_Tcell", 
                       "Bcell", 
                       "CD8_Tcell", 
                       "FCGR3A_monocyte", 
                       "NKcell", 
                       "DC", 
                       "platelet")

  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  pbmc$cell_type <- Idents(pbmc)

  sce <- as.SingleCellExperiment(pbmc)
  logcounts(sce) <- as.matrix(logcounts(sce))
  counts(sce) <- as.matrix(counts(sce))

} else {
  stop("Please pick either 'zeisel' or 'pbmc3k' for the dataset to base the simulation on.")
}


params <- splatEstimate(sce)

# deprob <- c(0.02, 0.06, 0.1)
# defacloc_defacscale <- c(0.75, 1.5)
# ncells = c(100, 1000, 10000, 100000, 20000, 50000, 80000, 1000000, 500000, 1000000, 500000)
ncells = c( 1000, 10000, 20000, 40000, 60000, 80000, 100000)
# ngenes = c(10000)
ngenes = c(2000)
# for (p in deprob){
for (ngene in ngenes){
  for (ncell in ncells){
    
    sim <- splatSimulate(params,
                         nGenes = ngene,
                         batchCells = ncell,
                         group.prob = c(0.25, 0.1, 0.1, 0.2, 0.3, 0.05),
                         method = "groups", 
                         de.prob = 0.1,
                         de.facLoc = 1.5,
                         de.facScale = 1.5,
                         out.prob = 0.001,
                         de.downProb = c(0),      
                         dropout.type = "experiment",
                         verbose = FALSE,
                         seed = 1234)
    
    
    clust.sim<- quickCluster(sim)
    sim <- computeSumFactors(sim, cluster=clust.sim, min.mean=0.1)
    sim <- logNormCounts(sim)
    sim <- runUMAP(sim)
    
    pumap <- plotUMAP(sim, colour_by = "Group")
    
    name <- paste0("ncell-", ncell,
                   "_ngene-", ngene)
    
    outdir <- file.path(maindir, name)
    dir.create(outdir)
    
    ggsave(plot = pumap, 
           filename =  paste0(outdir, "/", name, "_UMAP.png"))
    
    saveRDS(sim, paste0(outdir, 
                        "/",
                        name,
                        ".rds" ))
    
    png(paste0(outdir, "/", name, "_distribution.png"))
    grid = seq(0,30,.1)
    log_mean <- 1.5
    log_sd <- sqrt(1.5)
    
    plot(grid,dlnorm(grid, log_mean, log_sd),type="l",xlab="DE-factor",ylab="density")
    legend("topright", paste0("log-norm: ", log_mean, " mean, ", log_sd, " sd"),lty=1,col=1)
    dev.off()
      print("done!")
      print(ncell)
  }
    
}

