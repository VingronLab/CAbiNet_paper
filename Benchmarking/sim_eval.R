
gene_stats <- function(splatter_sim, grp2clust){
  
  groups <- sort(unique(colData(splatter_sim)$Group))
  
  rd <- rowData(splatter_sim)
  
  outdf <- data.frame("Group" = groups)
  outdf_list <- list()
  
  plots_list <- list()
  
  for (i in seq_len(length(groups))){
    
    grp <- as.character(groups[i])
    cluster <- grp2clust[[grp]]
    defacgrp <-paste0("DEFac",grp) 
    
    upreg <- rd$Gene[rd[,defacgrp] > 1]
    downreg <- rd$Gene[rd[,defacgrp] < 1]
    degenes <- na.omit(rd$Gene[rd$cluster == cluster])
    
    up <- length(intersect(degenes, upreg))/length(degenes)
    down <- length(intersect(degenes, downreg))/length(degenes)
    de <- length(intersect(degenes, c(upreg, downreg)))/length(degenes)
    
    totup <- length(intersect(degenes, upreg))/length(upreg)
    totdown <- length(intersect(degenes, downreg))/length(downreg)
    totde <- length(intersect(degenes, c(upreg, downreg)))/length(c(upreg, downreg))
    
    upcutoff <- get_defac_cutoff(splatter_sim = splatter_sim,
                                 group = grp, 
                                 cluster = cluster,
                                 cutoff = 0.1,
                                 direction = "up")
    
    downcutoff <- get_defac_cutoff(splatter_sim = splatter_sim,
                                   group = grp, 
                                   cluster = cluster,
                                   cutoff = 0.1,
                                   direction = "down")
    
    outdf_list[[i]] <- c("perc_up_of_clust"=up,
                         "perc_down_of_clust"=down,
                         "perc_de_of_clust"=de,
                         "perc_upclust_of_totup"=totup,
                         "perc_downclust_of_totdown" = totdown,
                         "perc_clust_of_totde" = totde,
                         "upper_quantile" = upcutoff,
                         "lower_quantile" = downcutoff)
    
    plots_list[[i]] <- list("de_genes" = plot_de_genes(splatter_sim, grp, cluster),
                            "all_genes" = plot_all_genes(splatter_sim, grp, cluster))
    
  }
  
  outdf <- data.frame("Group" = as.character(groups),
                      do.call(rbind, outdf_list))
  
  names(plots_list) <- groups
  
  return(list("stats" = outdf, "plots" = plots_list))
}


get_defac_cutoff <- function(splatter_sim, group, cluster, cutoff = 0.1, direction = "up"){
  
  defacgrp <-paste0("DEFac",group) 
  rd <- rowData(splatter_sim)
  
  if (direction == "up"){
    upreg_clust_idx <- which(rd[,defacgrp] > 1 & rd[,"cluster"] == cluster)
    factor_quant <- quantile(rd[upreg_clust_idx,defacgrp], cutoff)
    
  }
  
  if (direction == "down"){
    down_clust_idx <- which(rd[,defacgrp] < 1 & rd[,"cluster"] == cluster)
    factor_quant <- -quantile(rd[down_clust_idx, defacgrp], cutoff)
    
  }
  
  return(factor_quant)
  
}



get_defac_perc <- function(splatter_sim, group, cluster, direction = "up"){
  
  defacgrp <-paste0("DEFac",group) 
  rd <- rowData(splatter_sim)
  
  
  if (direction == "up"){
    upreg <- rd[,defacgrp][rd[,defacgrp] > 1]
    upreg_clust <- na.omit(rd[,defacgrp][rd[,defacgrp] > 1 & rd[,"cluster"] == cluster])
    upseq <- seq(1, max(upreg), by = 0.1)
    perc <- vector(mode="numeric", length = length(upseq))
    
    for (i in seq_along(upseq)){
      
      perc[i] <- sum(upreg_clust >= upseq[i])/sum(upreg >= upseq[i])
    }
    names(perc) <- upseq
  }
  
  if (direction == "down"){
    
    downreg <- rd[,defacgrp][rd[,defacgrp] < 1]
    downreg_clust <- na.omit(rd[,defacgrp][rd[,defacgrp] < 1 & rd[,"cluster"] == cluster])
    downseq <- seq(1, floor(min(downreg)), by = -0.1)
    
    perc <- vector(mode="numeric", length = length(downseq))
    
    for (i in seq_along(downseq)){
      
      perc[i] <- sum(downreg_clust >= downseq[i])/sum(downreg >= downseq[i])
    }
    names(perc) <- downseq
  }
  
  return(perc)
  
  
}


plot_de_genes <- function(splatter_sim, group, cluster){
  defacgrp <-paste0("DEFac",group) 
  rd <- rowData(splatter_sim)
  
  defac_df <- data.frame(
    "genes" = rd$Gene[rd[,defacgrp] > 1],
    "defac" = rd[,defacgrp][rd[,defacgrp] > 1],
    "clustered" = rd$Gene[rd[,defacgrp] > 1] %in% rd$Gene[rd$cluster == cluster]
  )
  
  ggplot(defac_df, aes(x=defac, fill = clustered)) +
    geom_histogram(alpha = 0.7) +
    theme_bw()
}

plot_all_genes <- function(splatter_sim, group, cluster){
  require(scales)
  
  defacgrp <-paste0("DEFac",group) 
  rd <- rowData(splatter_sim)
  
  defac_df_all <- data.frame("genes" = rd$Gene,
                             "defac" = rd[,defacgrp])
  defac_df_all$de <- "noDE"
  defac_df_all$de[rd[,defacgrp] > 1] <- "upregulated"
  defac_df_all$de[rd[,defacgrp] < 1] <- "downregulated"
  defac_df_all$cluster <- FALSE
  defac_df_all$cluster[rd$cluster == cluster] <- TRUE
  
  ggplot(defac_df_all, aes(x=defac, fill = cluster)) +
    geom_histogram(alpha = 0.7) +
    scale_y_continuous(trans = "pseudo_log") +
    theme_bw()
}



ari_cells <- function(reference, biclust_obj, reference_col = "Group"){
  
  cc <- biclust_obj@NumberxCol
  
  cc_clust <- vector(mode="numeric", length = ncol(cc))

  for (j in seq_len(ncol(cc))){

	  idx <- which(cc[,j] == TRUE)

	  if (length(idx) == 0) idx <- 0
#	  if (length(idx) > 1) idx <- idx[1]
	  if (length(idx) > 1) idx <- sample(x = idx, size = 1)
	  cc_clust[j] <- idx
  	
  }
  
  cc_clust <- as.factor(cc_clust)
  
  stopifnot(length(cc_clust) == length(colData(reference)[,reference_col]))
  
  ari <- mclust::adjustedRandIndex(cc_clust, colData(reference)[,reference_col])
  return(ari)
}

ari_cells_old <- function(reference, biclust_obj, reference_col = "Group"){
  
  cc <- biclust_obj@NumberxCol
  cc_clust <- as.numeric(seq_len(nrow(cc)) %*% cc)
  cc_clust <- as.factor(cc_clust)
  
  stopifnot(length(cc_clust) == length(colData(reference)[,reference_col]))
  
  ari <- mclust::adjustedRandIndex(cc_clust, colData(reference)[,reference_col])
  return(ari)
}

ari_genes <- function(splatter_sim, biclust_obj) {
  
  rd <- as.data.frame(rowData(splatter_sim))
  gc <- biclust_obj@RowxNumber
  
  rd <- rd[rd$Gene %in% rownames(gc),]
  
  de_fac <- dplyr::select(rd, starts_with("DEFacGroup"))
  de_fac <- as.matrix(de_fac)
  
  rMs <- rowMaxs(de_fac)
  idxs <- purrr::map(seq_len(nrow(de_fac)), function(x) which(de_fac[x,] == rMs[x]))
  
  de_genes <- de_fac > 1 | de_fac < 1 # doesnt consider that gene can be DE in 2 groups
  de_bool <- matrix(FALSE,
                    nrow = nrow(de_fac),
                    ncol = ncol(de_fac),
                    dimnames = list(rownames(de_fac),
                                    colnames(de_fac)))
  for (x in seq_len(nrow(de_fac))) {
    
    if(any(de_genes[x,])) {
    vals <- de_fac[x,][de_genes[x,]]
    max_val <- which(abs(vals) == max(abs(vals)))
    
    de_bool[x, which(de_fac[x,] == vals[max_val])] <- TRUE
    }else {
      next
    }
  }
  
  grp <- as.numeric(de_bool %*% seq_len(ncol(de_bool)))
  grp_fact <- as.factor(grp)
  names(grp_fact) <- rownames(de_bool)
  
  gc_clust <- as.numeric(gc %*% seq_len(ncol(gc)))
  gc_clust <- as.factor(gc_clust)
  names(gc_clust) <- rownames(gc)
  
  gc_clust <- gc_clust[order(base::match(names(gc_clust), names(grp_fact)))]
  
  ari <- mclust::adjustedRandIndex(gc_clust, grp_fact)
  
  return(ari)
}


logicM2vec = function(mat){
  # mat supposes to be logical matrix with clusters as rows, featurs as columns
  n = nrow(mat)
  vec = do.call(rbind, lapply(1:n, function(i){
    idx = which(mat[i,])
    df=data.frame(idx = idx)
    df$cluster = i
    return(df)}))
}




recovery.biclust = function(res, truth){
  temp = c()
  for (i in seq_len(truth@Number)){

    # Matrix genes x cells
    truthdf = Matrix::Matrix(0, 
                             nrow(truth@RowxNumber),
                             ncol(truth@NumberxCol), 
                             sparse = T)

    dimnames(truthdf) <- list(rownames(truth@RowxNumber), colnames(truth@NumberxCol))

    # Mark cells and genes in cluster i
    truthdf[truth@RowxNumber[,i], truth@NumberxCol[i,] ] <- 1
    
    if (res@Number > 1){
      
      # for each cluster in results, calculate U/I, take max for each cluster.
      temp <- c(temp, do.call(max, lapply(seq_len(res@Number), function(x){

        # Matrix marking cells & genes in results that are in cluster i.
        resdf <- Matrix::Matrix(0,
                               nrow(res@RowxNumber),
                               ncol(res@NumberxCol),
                               sparse = T)

        dimnames(resdf) <- list(rownames(res@RowxNumber),
                                colnames(res@NumberxCol))

        #Subset reference to rows & columns in results
        truth_tmp <- truthdf[rownames(truthdf) %in% rownames(resdf),
                             colnames(truthdf) %in% colnames(resdf)]

        resdf[res@RowxNumber[,x], res@NumberxCol[x,]] <- 1
        
        ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
        ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
        
        resdf <- resdf[ord_r, ord_c]
        
        # Intersection / Union
        return(sum(truth_tmp & resdf)/sum(truth_tmp|resdf))
      })))
      } else  {
        
        temp = c(temp, unlist(lapply(seq_len(res@Number), function(x){
          resdf = Matrix::Matrix(0, nrow(res@RowxNumber), ncol(res@NumberxCol), sparse = T)
          dimnames(resdf) <- list(rownames(res@RowxNumber), colnames(res@NumberxCol))
          resdf[res@RowxNumber[, x], res@NumberxCol[x,]] <- 1
          truth_tmp <- truthdf[rownames(truthdf) %in% rownames(resdf),
                               colnames(truthdf) %in% colnames(resdf)]
          
          ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
          ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
          
          resdf <- resdf[ord_r, ord_c]
          
          return(sum(truth_tmp & resdf)/sum(truth_tmp|resdf))
          
        })))
      }
  }
  return(sum(temp)/truth@Number)
}


relevance.biclust = function(res, truth){
  #column names of res/truth should be cluster, label
  temp = c()
  for (i in seq_len(res@Number)){
    
    resdf = Matrix::Matrix(0, nrow(res@RowxNumber), ncol(res@NumberxCol), sparse = T)
    dimnames(resdf) <- list(rownames(res@RowxNumber), colnames(res@NumberxCol))
    resdf[res@RowxNumber[,i], res@NumberxCol[i,]] <- 1
    
    if (res@Number > 0){
      
      temp = c(temp, do.call(max, lapply(seq_len(truth@Number), function(x){
        truthdf = Matrix::Matrix(0, nrow(truth@RowxNumber), ncol(truth@NumberxCol), sparse = T)
        dimnames(truthdf) <- list(rownames(truth@RowxNumber), colnames(truth@NumberxCol))
        truthdf[truth@RowxNumber[,x], truth@NumberxCol[x,] ] = 1
        truthdf <- truthdf[rownames(truthdf) %in% rownames(resdf),
                           colnames(truthdf) %in% colnames(resdf)]
        
        ord_r <- order(match(rownames(resdf), rownames(truthdf)))
        ord_c <- order(match(colnames(resdf), colnames(truthdf)))
        
        resdf <- resdf[ord_r, ord_c]
        
        return(sum(truthdf & resdf)/sum(truthdf|resdf))
        
      })))}else{
        
        temp = c(temp,  unlist(lapply(seq_len(truth@Number), function(x){
          truthdf = Matrix::Matrix(0, nrow(truth@RowxNumber), ncol(truth@NumberxCol), sparse = T)
          dimnames(truthdf) <- list(rownames(truth@RowxNumber), colnames(truth@NumberxCol))
          truthdf[truth@RowxNumber[,x], truth@NumberxCol[x,] ] = 1
          truthdf <- truthdf[rownames(truthdf) %in% rownames(resdf),
                             colnames(truthdf) %in% colnames(resdf)]
          
          ord_r <- order(match(rownames(resdf), rownames(truthdf)))
          ord_c <- order(match(colnames(resdf), colnames(truthdf)))
          
          resdf <- resdf[ord_r, ord_c]
          
          return(sum(truthdf & resdf)/sum(truthdf|resdf))
          
        })))
      }
  }
  return(sum(temp)/res@Number)
}



sim_truth <- function(splatter_sim,
                      factor_cutoff = 1,
                      no_overlap = TRUE){
  
  rd <- as.data.frame(rowData(splatter_sim))
  cd <- as.data.frame(colData(splatter_sim))
  
  de_fac <- dplyr::select(rd, starts_with("DEFacGroup"))
  colnames(de_fac) <- gsub("DEFacGroup", "Group", colnames(de_fac))
  
  de_fac <- as.matrix(de_fac)
  
  if (any(de_fac < 1)) stop("Only simulated datesets with no downregulated genes allowed.")
  de_genes <- de_fac > factor_cutoff
  
  if(isTRUE(no_overlap)){
    de_bool <- matrix(FALSE,
                      nrow = nrow(de_fac),
                      ncol = ncol(de_fac),
                      dimnames = list(rownames(de_fac),
                                      colnames(de_fac)))
    
    for (x in seq_len(nrow(de_fac))) {
      
      if(any(de_genes[x,])) {
        vals <- de_fac[x,][de_genes[x,]]
        max_val <- which(abs(vals) == max(abs(vals)))
        
        de_bool[x, which(de_fac[x,] == vals[max_val])] <- TRUE
      }else {
        next
      }
    }
    
    de_genes <- de_bool
  }
  
  cell_clusts <- model.matrix(~Group-1, data=cd) == 1
  colnames(cell_clusts) <- gsub("GroupGroup", "Group", colnames(cell_clusts))
  
  RowxNumber <- de_genes
  NumberxCol <- t(cell_clusts)
  
 if(isTRUE(no_overlap)){
   Number <- length(intersect(rownames(NumberxCol), colnames(RowxNumber)))
 } else {
   Number <- length(union(rownames(NumberxCol), colnames(RowxNumber)))
 }
 
  
  bic <- new("Biclust", "Parameters" = list("Splatter_Params" = splatter_sim@metadata$Params),
                       "RowxNumber" = RowxNumber,
                       "NumberxCol" = NumberxCol,
                       "Number" = Number,
                       "info" = list("splatter simulated data"))
  

}




evaluate_sim <- function(sce, biclust, truth_col = NULL){
  
  
  
  true_biclust <- sim_truth(splatter_sim = sce,
                            factor_cutoff = 1,
                            no_overlap = TRUE)
  
  if (biclust@Number > 0 ){
    
    nomono_biclust <- rm_monoclusters(biclust)
    nomono_truth <- rm_monoclusters(true_biclust)
    
    ac <- ari_cells(reference = sce,
                    biclust_obj = biclust,
                    reference_col = "Group")
    
	ac_old <- ari_cells_old(reference = sce,
							biclust_obj = biclust,
							reference_col = "Group")
    
    ag <- ari_genes(splatter_sim = sce,
                    biclust_obj = biclust)
    
    relevance <- relevance.biclust(nomono_biclust, nomono_truth)
    recovery <- recovery.biclust(nomono_biclust, nomono_truth)
	
	fARI <- fclust::ARI.F(VC = sce$Group, U = t(biclust@NumberxCol))

    genes_kept <- rownames(nomono_biclust@RowxNumber)
    nomono_truth@RowxNumber <- nomono_truth@RowxNumber[rownames(nomono_truth@RowxNumber) %in% genes_kept,]
    
    clustering_error <-  biclustlib_CE(nomono_biclust, nomono_truth)
    RNIA <-  biclustlib_RNIA(nomono_biclust, nomono_truth)

  }else{
    ac <- NA
    ag <- NA
    relevance = NA
    recovery = NA
    clustering_error = NA
    RNIA = NA
	fARI <- NA
	ac_old <- NA
  }
  
  evldf =  c("ARI_cells" = ac,
             "ARI_genes" = ag,
             "relevance" = relevance,
             "recovery" = recovery,
             "clustering_error" = clustering_error,
             "RNIA" = RNIA,
  			 "fuzzyARI_cells" = fARI,
  			 "ARI_cells_old" = ac_old)
  
  return(evldf)
  
}



evaluate_real <- function(sce, biclust, truth_col = NULL){
  
  if (biclust@Number > 0 ){

    nomono_biclust <- rm_monoclusters(biclust)

    ac <- ari_cells(reference = sce,
                    biclust_obj = biclust,
                    reference_col = truth_col)
    
    gini_res <- calgini(sce = sce,
                        biclust = biclust,
                        truth = truth_col,
                        name = "biclustering")


    ag <- NA
    relevance <- NA
    recovery <- NA
    relevance = NA
    recovery = NA
    clustering_error = NA
    RNIA = NA


  }else{
    ac <- NA
    ag <- NA
    relevance <- NA
    recovery <- NA
    clustering_error <- NA
    RNIA <- NA
    gini_res <- data.frame(detectg.mk = NA,
                           truthg.mk = NA,
                           silhouette = NA)

  }
  
  evldf <-  c("ARI_cells" = ac,
             "ARI_genes" = ag,
             "relevance" = relevance,
             "recovery" = recovery,
             "clustering_error" = clustering_error,
             "RNIA" = RNIA,
             "Gini_clusters" = gini_res$detectg.mk, # mean gini Idx per clustered genes
             "Gini_truth" = gini_res$truthg.mk,
             "Gene_silhoutte" = gini_res$silhouette,
              "Cell_silhouette" = gini_res$silhouette.cell
             )

  
  return(evldf)
  
}


name_biclust <- function(biclust, input){
  
  if (biclust@Number > 0){
    
    if (is.null(dimnames(biclust@RowxNumber))){
      dimnames(biclust@RowxNumber) <- list(rownames(input),
                                           paste0("BC", seq_len(ncol(biclust@RowxNumber))))
    }
    
    if (is.null(dimnames(biclust@NumberxCol))){

	  if (biclust@Number == 1 & 
        nrow(biclust@NumberxCol) != 1 & 
        ncol(biclust@NumberxCol) == 1){
        
   			biclust@NumberxCol <- t(biclust@NumberxCol)    
            
      }

      dimnames(biclust@NumberxCol) <- list(paste0("BC", seq_len(nrow(biclust@NumberxCol))),
                                           colnames(input))
    }
    
  
  }

  return(biclust)
}


qubic2biclust <- function(irisfgm_obj, sce, params){
  require(biclust)
  
  df_C = irisfgm_obj@BiCluster@CoCond_cell
  df_G = irisfgm_obj@BiCluster@CoReg_gene

  cc <- df_C$Condition
  names(cc) <- df_C$cell_name

  gc <- df_G$Condition
  names(gc) <- df_G$Gene

  ctypes = sort(unique(cc))
  gtypes = sort(unique(gc))
  bitypes = union(ctypes, gtypes)
  
  Number = length(bitypes)

  # stopifnot(identical(as.integer(bitypes), seq_len(Number)))

  # stopped here: You need to get the dimensions from the original sce and 
  # then only add the clusters on the ones actually coclustered.
  cells <- unique(names(cc))
  genes <- unique(names(gc))



  NumberxCol <- matrix(FALSE, 
                       ncol = ncol(sce),
                       nrow = Number)

  RowxNumber <- matrix(FALSE, 
                       ncol = Number,
                       nrow = nrow(sce))

  rownames(RowxNumber) <- rownames(sce)
  colnames(RowxNumber) <- paste0("BC", bitypes)

  rownames(NumberxCol) <- paste0("BC", bitypes)
  colnames(NumberxCol) <- colnames(sce)


  # NumberxCol <- matrix(FALSE, 
  #                      ncol = length(cells),
  #                      nrow = Number)

  # RowxNumber <- matrix(FALSE, 
  #                      ncol = Number,
  #                      nrow = length(genes))

  if (Number == 0){
    NumberxCol = matrix(0)
    RowxNumber = matrix(0)

    bic <- new("Biclust", "Parameters" = params,
                       "RowxNumber" = RowxNumber,
                       "NumberxCol" = NumberxCol,
                       "Number" = Number,
                       "info" = list("Results of QUBIC2"))
  
    return(bic)

  }else{

    for (x in seq_along(cells)){

        colidx <- which(colnames(NumberxCol) == cells[x])
        pick <- unique(cc[which(names(cc) == cells[x])])

        NumberxCol[pick, colidx] <- TRUE

    }


    for (y in seq_along(genes)){

        rowidx <- which(rownames(RowxNumber) == genes[y])
        pick <- unique(gc[which(names(gc) == genes[y])])

        RowxNumber[rowidx,pick] <- TRUE

    }
  }
  


  bic <- new("Biclust", "Parameters" = params,
                       "RowxNumber" = RowxNumber,
                       "NumberxCol" = NumberxCol,
                       "Number" = Number,
                       "info" = list("Results of QUBIC2"))
  
  return(bic)
}

get_qubic2_clusts <- function(file, params){

  tmp.block <- readLines(paste0(file, ".chars.blocks"))

  # Conditions
  keyword <- "Conds"
  tmp.bc <- grep(keyword, tmp.block, value = TRUE)
  tmp.cel.module <- vapply(strsplit(tmp.bc, ":", 2), "[", 2, FUN.VALUE = "character")
  CONDS <- as.character()  # store the conditions
  label_C <- as.numeric()  # store the occurence of one condistions

  for (j in seq_len(length(tmp.cel.module))) {
      BCcond <- unlist(strsplit(tmp.cel.module[j], split = " "))
      BCcond <- BCcond[BCcond != ""]  # exclude the blank string
      CONDS <- c(CONDS, BCcond)
      label_C <- c(label_C, rep(j, length(BCcond)))
  }
  df_C <- data.frame(cell_name = CONDS, Condition = label_C)
  df_C$cell_name <- as.character(df_C$cell_name)

  #Genes
  keyword = "Genes"
  tmp.bc <- grep(keyword, tmp.block, value = TRUE)
  tmp.cel.module <- vapply(strsplit(tmp.bc, ":", 2), "[", 2, FUN.VALUE = "character")
  GENES <- as.character()  # store the conditions
  label_G <- as.numeric()  # store the occurence of one condistions

  for (j in seq_len(length(tmp.cel.module))) {
      BCgene <- unlist(strsplit(tmp.cel.module[j], split = " "))
      BCgene <- BCgene[BCgene != ""]  # exclude the blank string
      GENES <- c(GENES, BCgene)
      label_G <- c(label_G, rep(j, length(BCgene)))
  }
  df_G <- data.frame(Gene = GENES, Condition = label_G)

  tmp.df_G <- df_G
  tmp.gene.list <- as.character(df_G$Gene)
  tmp.gene.list <- gsub("_[0-9]$", "", tmp.gene.list)
  df_G$Gene <- tmp.gene.list

  BCres <- qubic2biclust(df_C = df_C,
                         df_G = df_G,
                         params = params)

  return(BCres)
}




recovery_diff_size = function(res, truth){

  temp = c()
  for (i in seq_len(truth@Number)){

    # Matrix genes x cells
    truthdf = Matrix::Matrix(0, 
                             nrow(truth@RowxNumber),
                             ncol(truth@NumberxCol), 
                             sparse = T)

    dimnames(truthdf) <- list(rownames(truth@RowxNumber), colnames(truth@NumberxCol))

    # Mark cells and genes in cluster i
    truthdf[truth@RowxNumber[,i], truth@NumberxCol[i,] ] <- 1
    
    if (res@Number > 1){
      
      # for each cluster in results, calculate U/I, take max for each cluster.
      temp <- c(temp, do.call(max, lapply(seq_len(res@Number), function(x){

        # Matrix marking cells & genes in results that are in cluster i.
        resdf <- Matrix::Matrix(0,
                                nrow(res@RowxNumber),
                                ncol(res@NumberxCol),
                                sparse = TRUE)

        dimnames(resdf) <- list(rownames(res@RowxNumber),
                                colnames(res@NumberxCol))


        resdf[res@RowxNumber[,x], res@NumberxCol[x,]] <- 1


        if (!identical(dim(resdf), dim(truthdf))){
          in_both <- rownames(truthdf) %in% rownames(resdf)                          
          truth_tmp <- rbind(truthdf[in_both,], truthdf[!in_both,])

        } else {
          truth_tmp <- truthdf
        }
  

        ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
        ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
        resdf <- resdf[ord_r, ord_c]

        if (!identical(dim(resdf), dim(truthdf))){
          extra_genes <- Matrix::Matrix(0,
                                        nrow = sum(!in_both),
                                        ncol = ncol(resdf))

          resdf <- rbind(resdf, extra_genes)

        }


        # Intersection / Union
        return(sum(truth_tmp & resdf)/sum(truth_tmp | resdf))
      })))
      } else  {
        
        temp = c(temp, unlist(lapply(seq_len(res@Number), function(x){

          resdf = Matrix::Matrix(0,
                                 nrow(res@RowxNumber),
                                 ncol(res@NumberxCol),
                                 sparse = T)

          dimnames(resdf) <- list(rownames(res@RowxNumber),
                                  colnames(res@NumberxCol))
          
          resdf[res@RowxNumber[, x], res@NumberxCol[x,]] <- 1
          
          if (!identical(dim(resdf), dim(truthdf))){
            in_both <- rownames(truthdf) %in% rownames(resdf)                          
            truth_tmp <- rbind(truthdf[in_both,], truthdf[!in_both,])

          } else {
            truth_tmp <- truthdf
          }
    

          ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
          ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
          resdf <- resdf[ord_r, ord_c]

          if (!identical(dim(resdf), dim(truthdf))){
            extra_genes <- Matrix::Matrix(0,
                                          nrow = sum(!in_both),
                                          ncol = ncol(resdf))

            resdf <- rbind(resdf, extra_genes)

          }
          
          return(sum(truth_tmp & resdf)/sum(truth_tmp|resdf))
          
        })))
      }
  }
  return(sum(temp)/truth@Number)
}


relevance_diff_size = function(res, truth){
  #column names of res/truth should be cluster, label
  temp = c()

  for (i in seq_len(res@Number)){
    
    resdf = Matrix::Matrix(0,
                           nrow(res@RowxNumber),
                           ncol(res@NumberxCol),
                           sparse = T)

    dimnames(resdf) <- list(rownames(res@RowxNumber),
                            colnames(res@NumberxCol))
    
    resdf[res@RowxNumber[,i], res@NumberxCol[i,]] <- 1
    
    if (res@Number > 0){
      
      temp = c(temp, do.call(max, lapply(seq_len(truth@Number), function(x){

        truthdf = Matrix::Matrix(0,
                                 nrow(truth@RowxNumber),
                                 ncol(truth@NumberxCol),
                                 sparse = T)
        
        dimnames(truthdf) <- list(rownames(truth@RowxNumber),
                                  colnames(truth@NumberxCol))
        
        truthdf[truth@RowxNumber[,x], truth@NumberxCol[x,]] <- 1

        if (!identical(dim(resdf), dim(truthdf))){
          in_both <- rownames(truthdf) %in% rownames(resdf)                          
          truth_tmp <- rbind(truthdf[in_both,], truthdf[!in_both,])

        } else {
          truth_tmp <- truthdf
        }

        ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
        ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
        resdf <- resdf[ord_r, ord_c]

        if (!identical(dim(resdf), dim(truthdf))){
          extra_genes <- Matrix::Matrix(0,
                                        nrow = sum(!in_both),
                                        ncol = ncol(resdf))

          resdf <- rbind(resdf, extra_genes)

        }
          
        
        return(sum(truth_tmp & resdf)/sum(truth_tmp | resdf))
        
      })))}else{
        
        temp = c(temp,  unlist(lapply(seq_len(truth@Number), function(x){

          truthdf = Matrix::Matrix(0,
                                   nrow(truth@RowxNumber),
                                   ncol(truth@NumberxCol),
                                   sparse = T)

          dimnames(truthdf) <- list(rownames(truth@RowxNumber),
                                    colnames(truth@NumberxCol))

          truthdf[truth@RowxNumber[,x], truth@NumberxCol[x,]] = 1

          if (!identical(dim(resdf), dim(truthdf))){
            in_both <- rownames(truthdf) %in% rownames(resdf)                          
            truth_tmp <- rbind(truthdf[in_both,], truthdf[!in_both,])

          } else {
            truth_tmp <- truthdf
          }

          ord_r <- order(match(rownames(resdf), rownames(truth_tmp)))
          ord_c <- order(match(colnames(resdf), colnames(truth_tmp)))
          resdf <- resdf[ord_r, ord_c]

          if (!identical(dim(resdf), dim(truthdf))){
            extra_genes <- Matrix::Matrix(0,
                                          nrow = sum(!in_both),
                                          ncol = ncol(resdf))

            resdf <- rbind(resdf, extra_genes)

          }
          
          
          return(sum(truthdf & resdf)/sum(truthdf|resdf))
          
        })))
      }
  }
  return(sum(temp)/res@Number)
}


evaluate_sim_diff_size <- function(sce, biclust, truth_col = NULL){
  
  
  
  true_biclust <- sim_truth(splatter_sim = sce,
                            factor_cutoff = 1,
                            no_overlap = TRUE)
  
  if (biclust@Number > 0 ){
    
    nomono_biclust <- rm_monoclusters(biclust)
    nomono_truth <- rm_monoclusters(true_biclust)
    
    ac <- ari_cells(splatter_sim = sce,
                    biclust_obj = biclust)
    
    ag <- ari_genes(splatter_sim = sce,
                    biclust_obj = biclust)
    
    relevance <- relevance_diff_size(nomono_biclust, nomono_truth)
    recovery <- recovery_diff_size(nomono_biclust, nomono_truth)
    
    genes_kept <- rownames(nomono_biclust@RowxNumber)
    nomono_truth@RowxNumber <- nomono_truth@RowxNumber[rownames(nomono_truth@RowxNumber) %in% genes_kept,]
    
    clustering_error <-  biclustlib_CE(nomono_biclust, nomono_truth)
    RNIA <-  biclustlib_RNIA(nomono_biclust, nomono_truth)
  }else{
    ac <- NA
    ag <- NA
    relevance = NA
    recovery = NA
    clustering_error = NA
    RNIA = NA
  }
  
  evldf =  c("ARI_cells" = ac,
             "ARI_genes" = ag,
             "relevance" = relevance,
             "recovery" = recovery,
             "clustering_error" = clustering_error,
             "RNIA" = RNIA)
  
  return(evldf)
  
}


meanclust.sce <- function(sce, label){

    if (!(is.character(label))){
#         stop('label should be in colData of sce!')
        sce[['label']] = label
        label = 'label'
    }
    
    labels = sort(unique(sce[[label]]))
    df = matrix(0, nrow(sce), length(labels))

    for (i in 1:length(labels)){
        idx = which(sce[[label]] == labels[i])
        df[, i] = Matrix::rowMeans(counts(sce[, idx]))
    }

    colnames(df) = labels
    rownames(df) = rownames(sce)

    return(df)
}
meanclust.mat <- function(mat, label){
    ## label: a vector, clustering result of cells 
    labels = sort(unique(label))
    df = matrix(0, nrow(mat), length(labels))
    for (i in 1:length(labels)){
        idx = which(label == labels[i])
        df[, i] = Matrix::rowMeans(mat[, idx])
    }
    colnames(df) = labels
    rownames(df) = rownames(mat)
    return(df)
}


calgini = function(sce, biclust, truth = NULL, name = NULL){

    # cluster: clusters of cells
    cluster = rep(0, length = ncol(biclust@NumberxCol))

    for ( i in seq_len(nrow(biclust@NumberxCol))){
        cluster[biclust@NumberxCol[i,]] <- i
    }

    # mkgenes: genes in cluster
    mkgenes <- biclust@RowxNumber[rowSums(biclust@RowxNumber)>0,]

    ngenes <- nrow(mkgenes)
    nclust <- ncol(mkgenes)

    if(biclust@Number == 1){
      ngenes <- length(mkgenes)
      nclust <- 1
    }

    cluster.genes <- rep(0, length = ngenes)

    if(nclust == 1){
    
      cluster.genes[mkgenes] <- 1
    
    } else {

      for ( i in seq_len(nclust)){
        cluster.genes[mkgenes[,i]] = i
      }
    }



    if (is(sce, 'SingleCellExperiment')){
        
        # mean expression of genes per cluster
        detectm = meanclust.sce(sce, cluster)
        detectm = matrix(detectm, nrow = nrow(sce))
#         distmat = as.matrix(dist(counts(sce)))
        mat = counts(sce)
        
    }else if (is(sce, 'matrix')){
        mat = sce
        detectm = meanclust.mat(sce, cluster)
        if(is.null(dim(detectm))){
            stop('not applicable for a single cluster')
        }
#         distmat = as.matrix(dist(sce))
    }

    if (is.null(truth)){
        truthm = matrix(0)
        truthg.mk = matrix(0)
        truthg.else = matrix(0)

    } else {

        if (is(sce, 'SingleCellExperiment')){

          # mean expression of genes per real cluster
          truthm = meanclust.sce(sce, truth)
        } else if (is(sce, 'matrix')){
          truthm = meanclust.mat(sce, truth)
        }

        if (is.null(dim(truthm))){
            stop('not applicable for a single cluster')
        }

        # Gini index per clustered gene for true clusters
        idx = which(rowSums(biclust@RowxNumber) > 0)
        truthg.mk = apply(truthm[idx, ],  1, DescTools::Gini, na.rm = TRUE) 
        
        # Gini index for not clustered genes for true clusters
        idx = which(rowSums(biclust@RowxNumber) == 0)
        truthg.else = apply(truthm[idx, ],  1, DescTools::Gini , na.rm = TRUE)
    }
    
    # Gini index for clustered genes per cluster
    idx = which(rowSums(biclust@RowxNumber) > 0)

    # if(ncol(detectm) == 1){na.omit(DescTools::Gini(detectm[idx,]))}

    if (ncol(detectm) == 1) {
  
      detectg.mk <- purrr::map_dbl(detectm[idx, ], function(x) DescTools::Gini(x))

    } else {

      detectg.mk = apply(detectm[idx, ], 1, DescTools::Gini, na.rm = TRUE)

    }

    
    if (is(sce, 'SingleCellExperiment')){
        # distance between clustered genes (??)
        dist.mk = dist(counts(sce[idx,]))
    }else if (is(sce, 'matrix')){
        dist.mk = dist(sce[idx,])
    }
    
    # silhoutte score of clustered genes
    sil = cluster::silhouette(cluster.genes, dist.mk)
    
    dist.cell = dist(t(mat))                               
    sil.cell = cluster::silhouette(cluster, dist.cell)
    

    # Gini index for NOT clustered genes per cluster
    idx = which(rowSums(biclust@RowxNumber) == 0)

    if (ncol(detectm) == 1) {
  
      detectg.else <- purrr::map_dbl(detectm[idx, ], function(x) DescTools::Gini(x))

    } else {

      detectg.else = apply(matrix(detectm[idx, ]), 1, DescTools::Gini, na.rm = TRUE)

    }
    
    ginidf = data.frame(dataset = name,
                        detectg.mk = mean(detectg.mk), 
                        detectg.else = mean(detectg.else),
                        truthg.mk = mean(truthg.mk),
                        truthg.else = mean(truthg.else),
                        silhouette = if(all(is.na(sil))) NA else mean(sil[,3]),
                            silhouette.cell = if(all(is.na(sil.cell))) NA else mean(sil.cell[,3])
                       )
    
    # res = list(detectm = detectm,
    #            truthm = truthm, 
    #            detectg.mk = detectg.mk,
    #            detectg.else = detectg.else,
    #            truthg.mk = truthg.mk,
    #            truthg.else = truthg.else,
    #            ginidf = ginidf,
    #            silhouette = sil)

    return(ginidf)
    
}


de_comp <- function(biclust, ref_sce, truth_col, topdegs = 200){

  
      colLabels(ref_sce) <- factor(colData(ref_sce)[,truth_col])
      degs <- scoreMarkers(ref_sce)

      gcs <- biclust@RowxNumber
      gc_clusters <- colnames(gcs)

      intermat <- matrix(NA,
                         nrow = length(degs),
                         ncol = ncol(gcs))

      rownames(intermat) <- names(degs)
      colnames(intermat) <- colnames(gcs)


      for (ct in seq_len(length(degs))){

        degenes <- degs[[ct]]
        degenes <- degenes[order(degenes$min.logFC.cohen, decreasing = TRUE),]
        degenes <- degenes[seq_len(topdegs),]
        degenes <- as.data.frame(degenes) %>% 
                      rownames_to_column("gene")

        inter <- vector(mode="numeric", length = ncol(gcs))

        for (n in seq_len(ncol(gcs))){
            marker_genes <- rownames(gcs)[which(gcs[,n] == TRUE)]
            overlap_genes <- length(base::intersect(degenes$gene, marker_genes))
            inter[n] <- overlap_genes/topdegs
        }

        intermat[ct,] <- inter

      }

      return(intermat)

}


plot_de_intersection <- function(intermat){

  intermat <- as.data.frame(intermat)
  cnames <- colnames(intermat)
  p <- intermat %>% 
          rownames_to_column("cell_type") %>%
          pivot_longer(
              cnames,
              names_to = "cluster",
              values_to = "perc_intersection") %>%
          ggplot(aes(x=cluster, y = cell_type, fill = perc_intersection)) +
              geom_tile() +
              viridis::scale_fill_viridis(limits = c(0, 1))+
              theme_bw() + 
              theme(axis.text.x=element_text(angle=45, hjust=1))

  return(p)
}
