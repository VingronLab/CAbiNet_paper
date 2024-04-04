
source("./setup.R")

algorithm <- "CCA"


# CCA
cat("\nStarting CCA\n")

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
cat("\nFinished benchmarking!\n")
