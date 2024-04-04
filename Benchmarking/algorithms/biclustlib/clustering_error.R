library(reticulate)

source_python("./algorithms/biclustlib/subspace.py")


bic2list <- function(biclusters){
    bic_list <- list()
    for (i in seq_len(biclusters@Number)){
        bic_list[[i]] <- list(as.integer((which(biclusters@RowxNumber[,i]) - 1)),
                              as.integer((which(biclusters@NumberxCol[i,]) - 1)))
    }
    return(bic_list)
}


biclustlib_CE <- function(pred_biclusters, true_biclusters){
    
    nrows <- as.integer(nrow(pred_biclusters@RowxNumber))
    ncols <- as.integer(ncol(pred_biclusters@NumberxCol))
    
    ord_r <- order(match(rownames(pred_biclusters@RowxNumber),
                         rownames(true_biclusters@RowxNumber)))
    
    pred_biclusters@RowxNumber <- pred_biclusters@RowxNumber[ord_r, , drop = FALSE]

    
    ord_c <- order(match(colnames(pred_biclusters@NumberxCol),
                         colnames(true_biclusters@NumberxCol)))
    

    pred_biclusters@NumberxCol <- pred_biclusters@NumberxCol[, ord_c, drop = FALSE]
    
    

    pred_biclist <- bic2list(pred_biclusters)
    true_biclist <- bic2list(true_biclusters)
    
    pred_biclist <- convert2biclustlib(pred_biclist)
    true_biclist <- convert2biclustlib(true_biclist)
    
    ce <- clustering_error(pred_biclist, true_biclist, nrows, ncols)
    return(ce)
    
}

biclustlib_RNIA <- function(pred_biclusters, true_biclusters){
    
    nrows <- as.integer(nrow(pred_biclusters@RowxNumber))
    ncols <- as.integer(ncol(pred_biclusters@NumberxCol))
    

    ord_r <- order(match(rownames(pred_biclusters@RowxNumber),
                         rownames(true_biclusters@RowxNumber)))
    
    pred_biclusters@RowxNumber <- pred_biclusters@RowxNumber[ord_r, , drop = FALSE]
    
    ord_c <- order(match(colnames(pred_biclusters@NumberxCol),
                         colnames(true_biclusters@NumberxCol)))

    
    pred_biclusters@NumberxCol <- pred_biclusters@NumberxCol[, ord_c, drop = FALSE]
    
    pred_biclist <- bic2list(pred_biclusters)
    true_biclist <- bic2list(true_biclusters)
    
    pred_biclist <- convert2biclustlib(pred_biclist)
    true_biclist <- convert2biclustlib(true_biclist)
    
    RNIA <- relative_non_intersecting_area(pred_biclist, true_biclist, nrows, ncols)
    return(RNIA)
    
}



