

source("../setup.R")

algorithm <- "QUBIC2"


tmp_dir <- system("echo $MXQ_JOB_TMPDIR")
count_matrix <- file.path(tmp_dir,"count_matrix.tsv")
write.table(cnts, file = count_matrix, quote = FALSE, sep = '\t')

cmd <- paste("/home/kohl/bin/qubic -i", count_matrix, "-R -q", qQubic2, "-o", nclust, objF, sep=" ")


# QUBIC2
cat("\nStarting QUBIC2.\n")

t = Sys.time()

out <- system(cmd)

t.run = difftime(Sys.time(), t, units = 'secs')

res <- get_qubic2_clusts(file = count_matrix,
                         params = list("Call" = cmd,
                                       "-q" = qQubic2,
                                       "obj_fun" = objF))

##########
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
