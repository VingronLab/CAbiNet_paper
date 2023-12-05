
source("../setup.R")

algorithm <- "IRISFGM"

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


tmp_dir <- system("echo $MXQ_JOB_TMPDIR", intern = TRUE)

old_wd <- getwd()
setwd(tmp_dir)

iris_obj <- CreateIRISFGMObject(as.matrix(counts(data)))
iris_obj@Processed_count <- cnts



# QUBIC
cat("\nStarting IRIS-FGM.\n")

t = Sys.time()

iris_obj <- RunLTMG(iris_obj, Gene_use = "all", k = 5)
iris_obj <- CalBinaryMultiSignal(iris_obj)

iris_obj <- RunBicluster(iris_obj,
                    DiscretizationModel = "LTMG",
                    OpenDual = FALSE,
                    NumBlockOutput = nclust,
                    BlockOverlap = Qoverlap,
                    Extension = Qconsistency,
                    BlockCellMin = Qcmin)

t.run = difftime(Sys.time(), t, units = 'secs')

setwd(old_wd)

res <- qubic2biclust(irisfgm_obj=iris_obj,
                     sce = data,
                     params = opt)


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
