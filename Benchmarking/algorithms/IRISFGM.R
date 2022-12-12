
# WORK IN PROGRESS - DO NOT RUN!
# I did a test on this function, it works quite well. But I'm curious how the paramter 'k' in RunLTMG influences the results.
# How about varying the value of k too ?

source("/project/CAclust_scripts/CAclust_paper/benchmarking_splatter_sim/setup.R")

algorithm <- "IRISFGM"

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))
# data <- readRDS("/project/CAclust/data/tabula_muris.rds")

tmp_dir <- system("echo $MXQ_JOB_TMPDIR", intern = TRUE)

old_wd <- getwd()
setwd(tmp_dir)

iris_obj <- CreateIRISFGMObject(as.matrix(counts(data)))
iris_obj@Processed_count <- cnts
# iris_obj_p <- ProcessData(iris_obj, normalization = "cpm", IsImputation = FALSE)



# QUBIC
cat("\nStarting IRIS-FGM.\n")

t = Sys.time()

iris_obj <- RunLTMG(iris_obj, Gene_use = "all", k = 5)
iris_obj <- CalBinaryMultiSignal(iris_obj)
# iris_obj <- RunBicluster(iris_obj,
#                     DiscretizationModel = "LTMG",
#                     OpenDual = FALSE,
#                     NumBlockOutput = nclust)

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

# res <- name_biclust(biclust = res, input = cnts)

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
