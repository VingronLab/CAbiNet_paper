
source("../setup.R")



algorithm <- "Monocle3"



cds = new_cell_data_set(counts(data),
                         cell_metadata = colData(data))

# monocle
cat("\nStarting Monocle3 clustering\n")

t <- Sys.time()


cds <- preprocess_cds(cds,
                      num_dim = dims,
                      norm_method = "log",
                      use_genes = rownames(data))

cds <- reduce_dimension(cds,
                       reduction_method = reduction_method
                       )

cds <- cluster_cells(cds,
                     reduction_method = reduction_method,
                     cluster_method = 'leiden',
                     k = NNs,
                     resolution = resolution,
                     random_seed = 1)

part <- partitions(cds, reduction_method = reduction_method)



if (length(unique(part)) == 1){

    res <- biclust::BiclustResult(mypara = list(),
                              a = matrix(),
                              b = matrix(),
                              c = 0,
                              d = list())

    t.run <- difftime(Sys.time(), t, units = 'secs')

} else {

    marker_test_res <- top_markers(cds,
                                   group_cells_by = "partition",
                                   reduction_method = reduction_method,
                                   reference_cells = NULL,
                                   cores = 1,
                                   genes_to_test_per_group = genes_to_test_per_group,
                                   speedglm.maxiter = 200)

    top_specific_markers <- marker_test_res %>%
                                filter(fraction_expressing >= 0.10) %>%
                                filter(marker_test_q_value < 0.05) %>%
                                filter(marker_score > 0.1) %>%
                                group_by( gene_id) %>%
                                top_n(1, pseudo_R2)



    t.run <- difftime(Sys.time(), t, units = 'secs')

    #########

    ccs <- sort(unique(part))
    nclusts <- length(ccs)

    monocle_cells <- matrix(FALSE, nrow = length(ccs), ncol = ncol(cds))
    rownames(monocle_cells) <- paste0("Bic_", ccs)
    colnames(monocle_cells) <- Cells(cds)


    for (i in seq_along(ccs)){
        clust_cells <- Cells(cds)[which(part == ccs[i])]
        idx <- which(colnames(monocle_cells) %in% clust_cells)
        monocle_cells[i, idx] <- TRUE
    }

    gcs <- sort(unique(top_specific_markers$cell_group))

    monocle_genes <- matrix(FALSE, nrow = nrow(cds), ncol = length(ccs))
    rownames(monocle_genes) <- rownames(cds)
    colnames(monocle_genes) <- paste0("Bic_", ccs)

    for (i in seq_along(ccs)){

        if(!ccs[i] %in% gcs) next

        clust_genes <- top_specific_markers$gene_id[which(top_specific_markers$cell_group == ccs[i])]
        idx <- which(rownames(monocle_genes) %in% clust_genes)
        monocle_genes[idx, i] <- TRUE

    }

    res <- new("Biclust","Parameters" = opt,
                       "RowxNumber" = monocle_genes,
                       "NumberxCol" = monocle_cells,
                       "Number" = nclusts,
                       "info" = list("monocle clustering and DEA."))
}





#########
if (isTRUE(sim)){

    eval_res <-  evaluate_sim(sce = data_old,
                              biclust = res,
                              truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = nrow(cds),
                       "ncells" = ncol(cds),
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


} else{

    eval_res <-  evaluate_real(sce = data_old,
                            biclust = res,
                            truth_col = truth)

    eval_res <- c(list("algorithm" = algorithm),
                  as.list(eval_res),
                  list("ngenes" = nrow(cds),
                       "ncells" = ncol(cds),
                       "nclust_found" = res@Number,
                       "runtime" = t.run,
                       "runtime_dimreduc" = NA))

    eval_res <- bind_cols(eval_res, as_tibble(opt))


}


write_csv(eval_res, file.path(outdir, paste0(algorithm, "_", name, '_EVALUATION.csv')))

print('All done!')
# if (isTRUE(DECOMP)){

# intermat <- de_comp(biclust = res,
#                     ref_sce = data,
#                     truth_col = truth,
#                     topdegs = 100)

# p <- plot_de_intersection(intermat)

# ggsave(plot = p,
#        file = file.path(outdir, "img", paste0(algorithm, "_", name, '_DE_genes_comparison.png')))
# }

cat("\nFinished benchmarking!\n")
