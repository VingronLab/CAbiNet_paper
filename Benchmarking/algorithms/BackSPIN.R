
source("../setup.R")

# to install run "pip install backspinpy"
source("./backspin/backspin_fun.R")

algorithm <- "BackSPIN"


##################################
cat("\nStarting backSPIN.\n")
##################################


t = Sys.time()

#' @param data: 2-D array
#' the data matrix, rows should be genes and columns single cells/samples
#' @param numLevels: int
#' the number of splits that will be tried
#' @param first_run_iters: float
#' the iterations of the preparatory SPIN
#' @param first_run_step: float
#' the step parameter passed to _generate_widlist for the preparatory SPIN
#' @param runs_iters: int
#' the iterations parameter passed to the _divide_to_2and_resort.
#' influences all the SPIN iterations except the first
#' @param runs_step: float
#' the step parameter passed to the _divide_to_2and_resort.
#' influences all the SPIN iterations except the first
#' @param wid: float
#' the wid of every iteration of the splitting and resorting
#' @param split_limit_g: int
#' If the number of specific genes in a subgroup is smaller than this number
#' splitting of that subgrup is not allowed
#' @param split_limit_c: int
#' If the number cells in a subgroup is smaller than this number splitting of
#' that subgrup is not allowed
#' @param stop_const: float
#' minimum score that a breaking point has to reach to be suitable for splitting
#' @param low_thrs: float
#' genes with average lower than this threshold are assigned to either of the
#' splitting group reling on genes that are higly correlated with them
bs <- backspin(data = cnts,
               numLevels = numL, #2
               first_run_iters = 10.0,
               first_run_step = 0.05,
               runs_iters = 8,
               runs_step = 0.25,
               split_limit_g = 2,
               split_limit_c = 2,
               stop_const = stopc, #1.15
               low_thrs = lowT, # 0.2
               verbose = FALSE)

t.run = difftime(Sys.time(), t, units = 'secs')

#########

names(bs) <- c("genes_order", "cells_order", "genes_gr_level", "cells_gr_level", "cells_gr_level_sc", "genes_bor_level", "cells_bor_level") # nolint


# cell binary matrix
lvl <- numL + 1
cell_clusters <- bs$cells_gr_level[,lvl]
ccs <- unique(cell_clusters)

cell_mat <- matrix(FALSE, nrow = length(ccs), ncol = length(cell_clusters))
rownames(cell_mat) <- paste0("Bic_", ccs)
colnames(cell_mat) <- colnames(cnts)

for (i in seq_along(ccs)){
    idx <- bs$cells_order[which(cell_clusters == ccs[i])]  + 1
    cell_mat[i, idx] <- TRUE
}

# gene binary matrix

gene_clusters <- bs$genes_gr_level[,lvl]
gcs <- sort(unique(gene_clusters))

gene_mat <- matrix(FALSE, nrow = length(gene_clusters), ncol = length(gcs))
rownames(gene_mat) <- rownames(cnts)
colnames(gene_mat) <- paste0("Bic_", ccs)


for (i in seq_along(gcs)){
    idx <- bs$genes_order[which(gene_clusters == gcs[i])] + 1
    gene_mat[idx, i] <- TRUE
}


nclusts <- min(length(ccs), length(gcs))

# make biclust object

res <- new("Biclust","Parameters" = opt,
                       "RowxNumber" = gene_mat,
                       "NumberxCol" = cell_mat,
                       "Number" = nclusts,
                       "info" = list(paste0("BackSpin biclustering at level ", lvl )))


res <- CAbiNet::rm_monoclusters(res)



#########

if (isTRUE(sim)){

    eval_res <-  evaluate_sim(sce = data_old,
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
    
    eval_res <-  evaluate_real(sce = data_old,
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


cat("\nFinished benchmarking!\n")


###############
###### OLD ####
###############

# writes out file and reads it back in

# cef.dat = cnts

# # based on https://www.researchgate.net/publication/327258862_Supplementary_Information/data/5b849c604585151fd1370297/sdata2018160-s2.pdf
# tmp_dir = "/scratch/local/kohl/"
# output.cef = file.path(tmp_dir, paste0(algorithm, "_", name, '.cef'))
# clust_res.cef = file.path(tmp_dir, paste0(algorithm, "_", name, '_clusters.cef'))

# cef = cef.dat

# cef = rbind(gene="", data.frame(well="", cef))

# cef.head=c("CEF", "0", "1", "1", nrow(cef.dat), ncol(cef.dat), "0")

# write.table(matrix(cef.head, nrow=1),
#           output.cef,
#           sep="\t",
#           row.names=F,
#           col.names=F,
#           quote=F)

# write.table(cef,
#           output.cef,
#           sep="\t",
#           row.names=T,
#           col.names=NA,
#           quote=F,
#           append = T)


# backspin = paste0("backspin -i ", output.cef, " -o ", clust_res.cef, " &> ", clust_res.cef, ".log" )

# system(backspin)


# cells <- readr::read_delim(clust_res.cef,
#          n_max = 3,
#          skip = 1,
#          delim = "\t",
#          show_col_types = FALSE) %>%
#                select(well:last_col())


# genes <- readr::read_delim(clust_res.cef,                                                 
#                skip = 5,
#                delim = "\t",
#                show_col_types = FALSE) %>% 
#                 select(starts_with(c("gene", "Level")))
