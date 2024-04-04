

source("./setup.R")

algorithm <- "BC_Spectral"

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))

cnts <- as.matrix(counts(data))


# BC Spectral
cat("\nStarting BC Spectral\n")

t <- Sys.time()

res <- biclust::biclust(cnts,
                        method = BCSpectral(),
                        normalization = "log",
                        numberOfEigenvalues = numEig,
                        minr = minr_BCS,
                        minc = minc_BCS,
                        withinVar = withinVar)
                        # n_best = nbest,
                        # n_clusters = nclusters)

t.run <- difftime(Sys.time(), t, units = 'secs')

res <- name_biclust(biclust = res, input = cnts)

print(res)

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
cat("\nFinished benchmarking!\n")
