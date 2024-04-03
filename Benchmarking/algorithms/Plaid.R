
source("./setup.R")

algorithm <- "Plaid"


## PLAID
cat("\nStarting PLAID\n")

t = Sys.time()

res <- biclust::biclust(cnts,
                        method = BCPlaid(),
                        cluster = "b",
                        verbose = FALSE,
                        row.release = rrelease,
                        col.release = crelease)

t.run = difftime(Sys.time(), t, units = 'secs')

res <- name_biclust(biclust = res, input = cnts)

##################

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
