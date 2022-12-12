

source("/project/CAclust_scripts/CAclust_paper/benchmarking_splatter_sim/setup.R")

algorithm <- "scPNMF"

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))


# use logcounts.
# How many features? ALL!
ngenes = 100
dims_pca = 50

log_cnts <- as.matrix(logcounts(data_old))

res_pnmf <- scPNMF::PNMFfun(X = log_cnts,
                            K = 15,
                            method="EucDist",
                            tol=1e-4,
                            maxIter=1000,
                            verboseN = TRUE)

W <- res_pnmf$Weight
S <- res_pnmf$Score


W_select <- scPNMF::basisSelect(W = W,
                                S = S,
                                X = log_cnts,
                                toTest = TRUE,
                                toAnnotate = FALSE,
                                mc.cores = 1)


ig <- getInfoGene(W_select,
                  M = ngenes,
                  by_basis = T,
                  return_trunW = TRUE,
                  dim_use = NULL)


sel_genes <- unique(unlist(ig$InfoGene))

RowxNumber <- matrix(FALSE, nrow = length(sel_genes), ncol = length(ig$InfoGene))
rownames(RowxNumber) <- gtools::mixedsort(sel_genes)
colnames(RowxNumber) <- names(ig$InfoGene)

for (basis in seq_along(ig$InfoGene)){
    genes <- ig$InfoGene[[basis]]

    RowxNumber[genes, basis] <- TRUE
}

data <- data_old[sel_genes,]
data <- runPCA(data, ncomponents = dims_pca)
nn.clusters <- clusterCells(data, use.dimred="PCA")

# cnn.clusters <- clusterCells(data,
#                              use.dimred="PCA",
#                              BLUSPARAM=SNNGraphParam(k=10,
#                                                      type="rank",
#                                                      cluster.fun="walktrap"))

# scPNMF
cat("\nStarting scPNMF\n")

t <- Sys.time()

res <- biclust::biclust(cnts,
                        method = BCCC(),
                        delta = delta,
                        alpha=alpha,
                        number=nclust)

t.run <- difftime(Sys.time(), t, units = 'secs')

if (res@Number == 1 & ncol(res@RowxNumber) == ncol(res@NumberxCol)){
    res@NumberxCol <- t(res@NumberxCol)
}

res <- name_biclust(biclust = res, input = cnts)


#########
if (isTRUE(sim)){

    eval_res <-  evaluate_sim(sce = data,
                              biclust = res,
                              truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = nrow(cnts),
                       "ncells" = ncol(cnts),
                       "nclust_found" = res@Number,
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


} else{

    eval_res <-  evaluate_real(sce = data,
                            biclust = res,
                            truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = nrow(cnts),
                       "ncells" = ncol(cnts),
                       "nclust_found" = res@Number,
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


}


write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))

print('All done!')
