
source("/project/CAclust/paper-revision0/simulated/new_submission_scripts/setup.R")

algorithm <- "Bimax"

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))




cat("\nStarting Bimax\n")

t <- Sys.time()
res <-  biclust::biclust(cnts,
                        method = BCBimax(),
                        minr = minr,
                        minc = minc,
                        # maxc = maxc,
                        number=nclust) 

t.run <- difftime(Sys.time(), t, units = 'secs')

res <- name_biclust(biclust = res, input = cnts)



#########
if (isTRUE(sim)){

    eval_res <-  evaluate_sim(sce = data,
                              biclust = res,
                              truth_col = truth)

    eval_res <- c("algorithm" = algorithm,
                  eval_res,
                  "runtime" = t.run,
                  "runtime_dimreduc" = NA)

    eval_res <- bind_cols(as.list(eval_res), as_tibble(opt))

  
} else{
    
    eval_res <-  evaluate_real(sce = data,
                            biclust = res,
                            truth_col = truth)

    eval_res <- c("algorithm" = algorithm,
                  eval_res,
                  "runtime" = t.run,
                  "runtime_dimreduc" = NA)

    eval_res <- bind_cols(as.list(eval_res), as_tibble(opt))

  
}


write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))
