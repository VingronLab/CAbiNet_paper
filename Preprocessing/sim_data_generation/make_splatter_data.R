library(splatter)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(scran)

set.seed(12345)

# Uncomment for splatter data based on Zeisel Brain data.
# maindir <- "../sim_data/simulated_data/zeisel/"
# sce <- readRDS("../sim_data/data/preprocessed/zeisel.rds")

maindir <- "../sim_data/simulated_data/pbmc3k/"
sce <- readRDS("../sim_data/data/preprocessed/pbmc3k.rds")

params <- splatEstimate(sce)

deprob <- c(0.02, 0.06, 0.1)
defacloc_defacscale <- c(0.75, 1.5)

for (p in deprob){
  for (l in defacloc_defacscale){
    
    sim <- splatSimulate(params,
                         nGenes = 10000,
                         batchCells = 1000,
                         group.prob = c(0.25, 0.1, 0.1, 0.2, 0.3, 0.05),
                         method = "groups", 
                         de.prob = p,
                         de.facLoc = l,
                         de.facScale = l,
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
    
    name <- paste0("dePROB-", gsub("\\.","_", p),
                   "_defacLOC-",gsub("\\.","_", l),
                   "_defacSCALE-", gsub("\\.","_", s))
    
    outdir <- paste0(maindir, name)
    dir.create(outdir)
    
    ggsave(plot = pumap, 
           filename =  paste0(outdir, "/", name, "_UMAP.png"))
    
    saveRDS(sim, paste0(outdir, 
                        "/",
                        name,
                        ".rds" ))
    
    png(paste0(outdir, "/", name, "_distribution.png"))
    grid = seq(0,30,.1)
    log_mean <- l
    log_sd <- sqrt(s)
    
    plot(grid,dlnorm(grid, log_mean, log_sd),type="l",xlab="DE-factor",ylab="density")
    legend("topright", paste0("log-norm: ", log_mean, " mean, ", log_sd, " sd"),lty=1,col=1)
    dev.off()
  }
}

