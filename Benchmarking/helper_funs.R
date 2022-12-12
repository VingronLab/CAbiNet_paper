
##############################################
# This file contains helper functions for    #
# the evaluation of the biclustering results #
##############################################

ari_cells <- function(reference, biclust_obj, reference_col = "Group"){
  
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
    
    ag <- ari_genes(splatter_sim = sce,
                    biclust_obj = biclust)
    
    relevance <- relevance.biclust(nomono_biclust, nomono_truth)
    recovery <- recovery.biclust(nomono_biclust, nomono_truth)
    
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
        mat = counts(sce)
        
    }else if (is(sce, 'matrix')){
        mat = sce
        detectm = meanclust.mat(sce, cluster)
        if(is.null(dim(detectm))){
            stop('not applicable for a single cluster')
        }
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

    if (ncol(detectm) == 1) {
  
      detectg.mk <- purrr::map_dbl(detectm[idx, ], function(x) DescTools::Gini(x))

    } else {

      detectg.mk = apply(detectm[idx, ], 1, DescTools::Gini, na.rm = TRUE)

    }

    
    if (is(sce, 'SingleCellExperiment')){
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

###################################
### Additional helper functions ###
###################################

## remove monoclusters with less than a certain number (cutoff) of genes.
rm_raremonoclusters <- function(bic, cutoff = 10){
    gclust = gene_clusters(bic)
    sumg = summary(gclust)
    cclust = cell_clusters(bic)
    sumc = summary(cclust)
    common = intersect(unique(gclust), unique(cclust))

    realg = unique(names(sumg)[sumg > cutoff])
    common = union(common, realg)
    common = names(sort(sumc[names(sumc) %in% common], decreasing = T))
    
    idx_g = which(gclust %in% common)
    idx_c = which(cclust %in% common)
    
    bic@cell_clusters = bic@cell_clusters[idx_c]
    bic@cell_clusters = droplevels(bic@cell_clusters)
    bic@cell_clusters = factor(bic@cell_clusters, levels = common)
    
    bic@gene_clusters = bic@gene_clusters[idx_g]
    bic@gene_clusters = droplevels(bic@gene_clusters)
    bic@gene_clusters = factor(bic@gene_clusters, levels = common)
    
    bic = reindex(bic)
    
    stopifnot(!is.null(names(bic@cell_clusters)))
        stopifnot(!is.null(names(bic@gene_clusters)))
    
     if(!is.empty(bic@SNN)){
          
          selr <- which(rownames(bic@SNN) %in% c(names(bic@cell_clusters),
                                                 names(bic@gene_clusters)))
          
          selc <- which(colnames(bic@SNN) %in% c(names(bic@cell_clusters),
                                                 names(bic@gene_clusters)))
          bic@SNN <- bic@SNN[selr, selc]
          
        }
    return(bic)
    
}

reindex = function(bic){
    clusters = unique(levels(gene_clusters(bic)))
    names(clusters) = seq(clusters) 

    newbic = bic
    newindex = clusters[gene_clusters(bic)]
    newbic@gene_clusters = factor(names(newindex), levels = seq(clusters))
    names(newbic@gene_clusters) = names(bic@gene_clusters)

    newindex = clusters[cell_clusters(bic)]
    newbic@cell_clusters = factor(names(newindex), levels = seq(clusters))
    names(newbic@cell_clusters) = names(bic@cell_clusters)

    return(newbic)
}

rm_raremonoclusters_sce <- function(sce, cutoff = 10){
    bic = metadata(sce)$caclust
    bic = rm_raremonoclusters(bic, cutoff = cutoff)
    metadata(sce)$caclust = bic
    sce = sce[ , colnames(sce) %in% names(bic@cell_clusters)]
    colData(sce)$caclust = cell_clusters(bic)
    return(sce)
}

features_biMAP <- function(sce,
                         caclust, 
                         features, 
                         color_cells_by="expression", 
                         assay = "logcounts",
                         gene_alpha = 0.5,
                         cell_alpha = 0.8,
                         gene_size = 1,
                         cell_size = 1,
                          min_cutoff = 0,
                          max_cutoff = NA,
                          nrow = 2,
                           ncol = NA,
                          show_legend = TRUE,
                           show_labs = TRUE,
                          filename = NULL,
                          width = 8,
                          height = 4){
  
  stopifnot(is(caclust, "caclust"))
  umap_coords <- caclust@bimap
  variables = colnames(umap_coords)
  
  stopifnot(is(sce, "SingleCellExperiment"))
  stopifnot(length(features)>=1)
  
  if(color_cells_by == "expression") {
    discr <- FALSE
    if(!is.null(features)) lgnd <- 'Expression'
  }else{
    discr <- TRUE
    lgnd <- color_cells_by
  }
  
  cell_idx <- which(umap_coords$type == "cell")
  
  if(sum(features %in% rownames(umap_coords)) < length(features) ){
    if(sum(features %in% rownames(umap_coords)) == 0){
      stop("All Features Not Found!")
    }
    features = features[features %in% rownames(umap_coords)]
    warning("Not all of the features are available.")
  }
  # stopifnot(isTRUE(feature %in% umap_coords$name))
  cnts <- SummarizedExperiment::assay(sce, assay)
  FetchData <- function(cnts, umap_coords, features){
      for (feature in features){
          umap_coords[ ,feature] <- NA
          umap_coords[cell_idx ,feature] <- cnts[feature, umap_coords$name[cell_idx]]
      }
      return(umap_coords)
  }
    
  umap_coords = FetchData(cnts = cnts, 
                          umap_coords = umap_coords, 
                          features = features)
    
    
   min_cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(umap_coords[cell_idx, feature]),
        no = cutoff
      ))
    },
    cutoff = min_cutoff,
    feature = features
  )
  max_cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(umap_coords[cell_idx, feature]),
        no = cutoff
      ))
    },
    cutoff = max_cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min_cutoff, max_cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
    
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  
  # Apply cutoffs
    if (length(features) > 1){
      umap_coords[, features] <- sapply(
      X =seq(length(features)), #which(featues %in% colnames(umap_coords))
      FUN = function(index) {
      data.feature <- as.vector(x =umap_coords[, features][,index])
          min.use = quantile(min_cutoff, 0.01)
          max.use = quantile(max_cutoff, 0.99)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    }
  )}
    
     # Make list of plots
  plots <- list()
    
    for(i in 1:length(features)){
        feature = features[i]
   umap_coords[cell_idx, ] <- umap_coords[cell_idx, ][order(umap_coords[cell_idx, feature], decreasing = FALSE),]
        
      plot = suppressWarnings({
      ggplot2::ggplot()+
        ggplot2::geom_point(umap_coords[umap_coords$type == "gene", ],
                   mapping=ggplot2::aes_(~x, ~y, text = paste0(
                                           "Type: ", quote(type), "\n",
                                           "Name: ", quote(name), "\n",
                                           "Cluster: ", quote(cluster))), 
                   color ="#A9A9A9", 
                   alpha = gene_alpha,
                   size = gene_size) +  #grey
        ggplot2::geom_point(umap_coords[umap_coords$type == "cell", ],
                   mapping=ggplot2::aes_(~x, 
                                         ~y, 
                                         color = as.name(feature), 
                                         text = paste0(
                                            "Type: ", quote(type), "\n",
                                            "Name: ", quote(name), "\n",
                                            "Cluster: ", quote(cluster))),
                   alpha = cell_alpha,
                   size = cell_size) + 
        ggplot2::geom_point(data = na.omit(umap_coords[feature,c("name", "x","y")]),
                            ggplot2::aes_(~x, ~y),
                   color = "red") +
        ggrepel::geom_text_repel(data = na.omit(umap_coords[feature,c("name", "x","y")]),
                                 ggplot2::aes_(~x, ~y, label= ~name),
                                 color = "red") +
        # ggplot2::scale_color_viridis_c(name=lgnd, discrete = discr, limits = c(0,9)) +
        ggplot2::scale_color_viridis_c(name=lgnd,  limits = c(floor(min(min_cutoff)), ceiling(max(max_cutoff)))) +
          
        ggplot2::labs(x="Dim 1",
                      y="Dim 2")+
        ggplot2::theme_bw() 
      })
        if (isFALSE(show_labs)){
            plot = plot + theme(axis.text.x = element_blank(), 
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank()
                               )
        }
        if (isTRUE(show_legend)){
            plots[[i]] = plot #+ theme(legend.position = 'none')
            }else{
            plots[[i]] = plot + theme(legend.position = 'none')
        }
        }

    p = grid.arrange( grobs = plots,  nrow = nrow, ncol = ncol)
    if (!is.null(filename)){
        ggsave(paste0(filename, '.pdf'), p, width = width , height = height)
        }
    return (p)
}
