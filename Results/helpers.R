suppressPackageStartupMessages({
  
    library(ComplexHeatmap)
    library(SingleCellExperiment)
    library(scater)
    library(scuttle)
    library(scran)
    library(RColorBrewer)
    library(tidyverse)
    library(ggplot2)
    library(mclust)
    library(ggrepel)
    library(gridExtra)
    library(cowplot)
    library(infotheo)
    library(stringr)
    library(cluster)
    library(entropy)
    library(sva)
    library(Seurat)
    library(umap)
    library(reticulate)
    library(CAbiNet)

})

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
    
     if(!CAbiNet:::is.empty(bic@SNN)){
          
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

## plot several feature_biMAP together and make the range of color codes in each feature_bimap in the same scale
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
  # brewer.gran <- ifelse(
  #   test = length(x = cols) == 1,
  #   yes = brewer.pal.info[cols, ]$maxcolors,
  #   no = length(x = cols)
  # )
  # umap_coords[cell_idx,] <- umap_coords[cell_idx,][order(umap_coords[cell_idx,"expression"], decreasing = FALSE),]
  
  # Apply cutoffs
    if (length(features) > 1){
      umap_coords[, features] <- sapply(
      X =seq(length(features)), #which(featues %in% colnames(umap_coords))
      FUN = function(index) {
      data.feature <- as.vector(x =umap_coords[, features][,index])
      # min.use <- SetQuantile(cutoff = min_cutoff[index], data.feature)
      # max.use <- SetQuantile(cutoff = max_cutoff[index], data.feature)
          min.use = quantile(min_cutoff, 0.01)
          max.use = quantile(max_cutoff, 0.99)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      # if (brewer.gran == 2) {
      #   return(data.feature)
      # }
      # data.cut <- if (all(data.feature == 0)) {
      #   0
      # }
      # else {
      #   as.numeric(x = as.factor(x = cut(
      #     x = as.numeric(x = data.feature),
      #     breaks = brewer.gran
      #   )))
      # }
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


## Heatmap function for visualizing biclusters
library(reticulate)
library(CAbiNet)
library(tidyverse)

library(SingleCellExperiment)
library(scater)
library(scran)

library(ComplexHeatmap)
library(viridis)
library(Matrix)
library(circlize)
library(biclust)
library(ggthemes)

library(s4vd)

source("..//Benchmarking/sim_eval.R")


get_gene_clusters <- function(splatter_sim) {

  rd <- as.data.frame(rowData(splatter_sim))

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

  return(grp_fact)
}

make_bic_hm <- function(bic, sim, max_bics = bic@Number, algorithm = "Replace with algorithm name"){
    
    nbics <- bic@Number
    max_bics <- min(max_bics, bic@Number)
    
    clust_order_cells <- c()
    for (b in seq_len(nbics)){
        idx <- which(bic@NumberxCol[b,] == TRUE)
        idx <- idx[!idx %in% clust_order_cells]
        clust_order_cells <- c(clust_order_cells, idx)    
    }
    
    idx <- seq_len(ncol(bic@NumberxCol))
    idx <- idx[!idx %in% clust_order_cells]
    clust_order_cells <- c(clust_order_cells, idx)
    
    
    clust_order_genes <- c()
    for (b in seq_len(nbics)){
        idx <- which(bic@RowxNumber[,b] == TRUE)
        idx <- idx[!idx %in% clust_order_genes]
        clust_order_genes <- c(clust_order_genes, idx)    
    }
    idx <- seq_len(nrow(bic@RowxNumber))
    idx <- idx[!idx %in% clust_order_genes]
    clust_order_genes <- c(clust_order_genes, idx)
    
    
    data_sort <- sim[clust_order_genes, clust_order_cells]
    gcs <- get_gene_clusters(splatter_sim = data_sort)
    gcs <- paste0("Group", gcs)
    
    pal <- tableau_color_pal("Tableau 10")
    annocol <- pal(nbics)
    
    col_df <- bic@NumberxCol
    col_df[bic@NumberxCol] <- "y"
    col_df[!bic@NumberxCol] <- "n"
    col_df <- as.data.frame(t(col_df))
    col_df <- col_df[clust_order_cells,, drop = FALSE]
#    col_df <- as.data.frame(col_df) # necessaey when max_bics==1
    col_df <- col_df[,seq_len(max_bics), drop=FALSE]
    
    col_df$Truth <- data_sort$Group
    col_df <- rev(col_df)
    
    anno_list <- list()
    for (b in seq_len(nbics)){
        anno_list[[paste0("BC",b)]] <- c("y" = annocol[b], "n" = "white")
    }
    anno_list <- anno_list[seq_len(max_bics)]
    
    truecol <- pal(length(unique(data_sort$Group)))
    names(truecol) <- sort(unique(data_sort$Group))
    
    anno_list_cols <- anno_list
    anno_list_cols[["Truth"]] <- truecol
    
    # bind_rows(anno_df)
    
    column_ha = HeatmapAnnotation(df = col_df,
                                  col = anno_list_cols,
                                  show_legend = c(FALSE), 
                                  # annotation_name_gp = gpar(fontsize = 3),
                                  na_col = "white")
    
    row_df <- bic@RowxNumber
    print(dim(row_df))
    row_df[bic@RowxNumber] <- "y"
    row_df[!bic@RowxNumber] <- "n"
    row_df <- as.data.frame(row_df)
    print(dim(row_df))
    row_df <- row_df[clust_order_genes,, drop = FALSE]
#	row_df <- as.data.frame(row_df)
    row_df <- row_df[,seq_len(max_bics), drop=FALSE]
    
    row_df$Truth <- gcs
    row_df <- rev(row_df)
    
    anno_list_rows <- anno_list
    anno_list_rows[["Truth"]] <- c(truecol, "Group0" = "white") 
    
    row_ha = rowAnnotation(df = row_df,
                           col = anno_list_rows,
                           show_legend = c(FALSE), 
                           annotation_name_gp = gpar(fontsize = 0),
                           na_col = "white")
    
    
    cnts <- as.matrix(logcounts(data_sort))
    
    col_fun = colorRamp2(seq(quantile(cnts, 0.01), quantile(cnts, 0.99), length = 10), viridis(10))
    
    hm <- Heatmap(matrix = cnts,
                      col = col_fun,
                      name = " ",
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      show_row_names = FALSE,
                      show_column_names = FALSE,
                      column_title = algorithm,
                      top_annotation = column_ha,
                      left_annotation = row_ha,
                      use_raster = TRUE)

    return(hm)
    
}
