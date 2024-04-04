## R wrapper for divbiclust cpp functions
source("./setup.R")
sourceCpp('./algorithms/src/divbiclust.cpp')

algorithm = 'divbiclust'

cat("\nStarting divbiclust...\n")
t = Sys.time()


res = DivBiclust(ds_type = ds_type, 
				 in_file = in_file,
          		 max_diff = max_diff, 
          		 do_rate = do_rate, 
          		 seed_col_sz = seed_col_sz, 
          		 max_col_sz = max_col_sz, 
          		 simthresh = simThresh)


t.run = Sys.time() - t

ari = stringr::word(res, 1, 1, '_')
ari = as.numeric(ari)
ncluster = stringr::word(res, 2, 2, '_')
ncluster = as.numeric(ncluster)

eval_res <- c(list("algorithm" = algorithm),
              'ARI_cells' = ari,
              # 'ARI_genes' = NA,
              # 'relevance' = NA,
              # 'recovery' = NA,
              # 'clustering_error' = NA,
              # 'RNIA' = NA,
              list("ngenes" = ntop,
                   "ncells" = ncell,
                   "nclust_found" = ncluster,
                   "runtime" = t.run,
                   "runtime_dimreduc" = NA))

eval_res <- bind_cols(eval_res, as_tibble(opt))



write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))
cat("\nFinished benchmarking!\n")
