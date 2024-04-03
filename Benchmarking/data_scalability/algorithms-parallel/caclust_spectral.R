
#source("/project/CAclust_scripts/CAclust_paper/benchmarking_splatter_sim/newARI/setup.R")
#source('../../setup.R')
source('/project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/Benchmarking/setup.R')
algorithm <- "caclust_spectral"

# write_csv(as_tibble(opt), file.path(outdir, paste0(algorithm, "_", name, '_parameters.csv')))
# saveRDS(opt, file.path(outdir, paste0(algorithm, "_", name, '_parameters.rds')))




# if (isTRUE(graph_select)){
#     ngenes <- nrow(counts) * 0.8
# }

cat("\nStarting CA.\n")
t = Sys.time()
caobj = cacomp(cnts,
               dims = dims,
               top = nrow(cnts),
               python = TRUE)

t.CA = difftime(Sys.time(), t, units = 'secs')


# Spectral clustering
cat("\nStarting CAclust spectral.\n")

t = Sys.time()

res <-
# tryCatch({

    caclust(obj = caobj,
                k = NNs,
                loops = FALSE,
                SNN_prune = prune,
                mode = SNN_mode,
                select_genes = graph_select,
                # prune_overlap = TRUE,
                prune_overlap = prune_overlap,
                overlap = overlap,
                calc_gene_cell_kNN = gcKNN,
                algorithm = 'spectral',
                spectral_method = "skmeans",
                use_gap = usegap,
                nclust = nclust,
                BPPARAM = MulticoreParam(workers = 16))
# },
# error = function(cond){
#
#     return(biclust::BiclustResult(mypara = list(),
#                                   a = matrix(),
#                                   b = matrix(),
#                                   c = 0,
#                                   d = list()))
#
# })

t.run = difftime(Sys.time(), t, units = 'secs')
# write.csv(t.SC, file.path(outdir, paste0(name, '_sc_runtime.csv')), row.names =F )
if (is(res, "caclust")) res <- convert_to_biclust(res)

if (isTRUE(sim)){


    eval_res <- c(list("algorithm" = algorithm),
                  list("ngenes" = nrow(cnts),
                       "ncells" = ncol(cnts),
                       "runtime" = t.run,
                       "runtime_dimreduc" = t.CA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))

 
} else{
   

    eval_res <- c(list("algorithm" = algorithm),
                  list("ngenes" = nrow(cnts),
                       "ncells" = ncol(cnts),
                       "runtime" = t.run,
                       "runtime_dimreduc" = t.CA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


}

write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))
print('All done!')
cat("\nFinished benchmarking!\n")
