

source("/project/CAclust_scripts/CAclust_paper/benchmarking_splatter_sim/setup.R")

algorithm <- "sv4d"

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))


# sv4d
cat("\nStarting sv4d\n")

t <- Sys.time()

res <- biclust::biclust(cnts,
                        method=BCs4vd(),
                        pcerv=pcerv,
                        pceru=pceru,
                        ss.thr=ss_thr,
                        start.iter=3,
                        size=0.632,
                        cols.nc=TRUE,
                        steps=100,
                        pointwise=TRUE,
                        merr=0.0001,
                        iter=100,
                        nbiclust=nclust,
                        col.overlap=FALSE)

t.run <- difftime(Sys.time(), t, units = 'secs')
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
