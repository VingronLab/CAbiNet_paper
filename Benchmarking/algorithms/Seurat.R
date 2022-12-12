
source("/project/CAclust_scripts/CAclust_paper/benchmarking_splatter_sim/setup.R")

algorithm <- "Seurat"

rownames(data_old) = gsub('_', '-', rownames(data_old))

# if (is(counts(data_old), 'dgTMatrix')){
#   counts(data_old) = Matrix(counts(data_old), sparse = T)
# }
# if (is(logcounts(data_old), 'dgTMatrix')){
#   logcounts(data_old) = Matrix(logcounts(data_old), sparse = T)
# }


# seu <- as.Seurat(data_old)

seu <- CreateSeuratObject(counts = as(counts(data_old), "dgCMatrix"),
                          meta.data = as.data.frame(colData(data_old)))


seu <- SetAssayData(object = seu,
                    slot = "data",
                    new.data = as(logcounts(data_old), "dgCMatrix"))
                    
seu <- FindVariableFeatures(object = seu,
                            nfeatures = ntop)
# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 4000)

all.genes <- rownames(seu)


# Seurat
cat("\nStarting Seurat clustering\n")

t <- Sys.time()


seu <- ScaleData(object = seu, features = VariableFeatures(seu))

seu <- RunPCA(object = seu,
              npcs = dims,
              features = VariableFeatures(seu))

seu <- FindNeighbors(object = seu,
                     dims = seq_len(dims),
                     k.param = NNs,
                     prune.SNN = prune)

seu <- FindClusters(object = seu,
                    resolution = resol)

seu.markers <- FindAllMarkers(seu,
                              features = VariableFeatures(seu),
                              only.pos = TRUE,
                              min.pct = min_perc, # 0.25
                              logfc.threshold = logfc_thr, # 0.25
                              return.thresh = rthr) # 0.01


t.run <- difftime(Sys.time(), t, units = 'secs')

#########

ccs <- sort(unique(seu$seurat_clusters))
nclusts <- length(ccs)

seurat_cells <- matrix(FALSE, nrow = length(ccs), ncol = ncol(seu))
rownames(seurat_cells) <- paste0("Bic_", ccs)
colnames(seurat_cells) <- Cells(seu)


for (i in seq_along(ccs)){
    clust_cells <- Cells(seu)[which(seu$seurat_clusters == ccs[i])]
    idx <- which(colnames(seurat_cells) %in% clust_cells)
    seurat_cells[i, idx] <- TRUE
}

gcs <- sort(unique(seu.markers$cluster))

seurat_genes <- matrix(FALSE, nrow = length(VariableFeatures(seu)), ncol = length(ccs))
rownames(seurat_genes) <- VariableFeatures(seu)
colnames(seurat_genes) <- paste0("Bic_", ccs)

# seurat_genes <- matrix(FALSE, nrow = nrow(seu), ncol = length(ccs))
# rownames(seurat_genes) <- rownames(seu)

for (i in seq_along(ccs)){

    if(!ccs[i] %in% gcs) next

    clust_genes <- seu.markers$gene[which(seu.markers$cluster == ccs[i])]
    idx <- which(rownames(seurat_genes) %in% clust_genes)
    seurat_genes[idx, i] <- TRUE

}



res <- new("Biclust","Parameters" = opt,
                       "RowxNumber" = seurat_genes,
                       "NumberxCol" = seurat_cells,
                       "Number" = nclusts,
                       "info" = list("Seurat clustering and DEA."))

data_old <- data_old[VariableFeatures(seu),]
#cat('length of variable features,', length(VariableFeatures(seu)), '\n')

#data_old <- data_old[rownames(data_old) %in% VariableFeatures(seu),]

#cat('Number of genes remained:', nrow(data_old), '\n')

#########
if (isTRUE(sim)){

    eval_res <-  evaluate_sim(sce = data_old,
                              biclust = res,
                              truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = length(VariableFeatures(seu)),
                       "ncells" = ncol(seu),
                       "nclust_found" = res@Number,
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


} else{

    eval_res <-  evaluate_real(sce = data_old,
                            biclust = res,
                            truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = length(VariableFeatures(seu)),
                       "ncells" = ncol(seu),
                       "nclust_found" = res@Number,
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


}


write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))

print('All done!')
