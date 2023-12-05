
library(Matrix)
library(Rcpp)
library(CAbiNet)
library(APL)
library(tidyverse)
library(SingleCellExperiment)
library(optparse)
library(reticulate)
library(scran)
library(scater)
library(IRISFGM)
library(Seurat)
library(mclust)
library(isa2)
library(s4vd)
library(QUBIC)
library(runibic)
library(skmeans)
library(scPNMF)
library(bluster)
# library(monocle3)

set.seed(2358)

# source("./helper_funs.R")
# source("./clustering_error.R")


option_list = list(
    make_option(c("--name"),
                type="character",
                action = "store",
                default="sample",
                help="Name of sample",
                metavar="character"),

    make_option(c("--file"),
                type="character",
                action = "store",
                default=NULL,
                help="Name of file to load",
                metavar="character"),

    make_option(c("--dataset"),
                type="character",
                action = "store",
                default=NULL,
                help="Name of the dataset",
                metavar="character"),

    make_option(c("--outdir"),
                type="character",
                action = "store",
                default=NULL,
                help="output directory",
                metavar="character"),

    make_option(c("--decomp"),
                type="logical",
                action = "store_true",
                default=FALSE,
                help="Whether DECOMP plot should be made",
                metavar="logical"),

    make_option(c("--ntop"),
                type="numeric",
                action = "store",
                default=NA,
                help="top X most variable genes",
                metavar="numeric"),

    make_option(c("--dims"),
                type="numeric",
                action = "store",
                default=NA,
                help="dimensions",
                metavar="numeric"),

    make_option(c("--vst"),
                type="logical",
                action = "store_true",
                default=FALSE,
                help="Whether to use vst or not",
                metavar="logical"),

    make_option(c("--distance"),
                type="character",
                action = "store",
                default='Euclidean',
                help="distance measure for gene-gene/sample-sample adjacency matrix",
                metavar="character"),

    make_option(c("--NNs"),
                type="numeric",
                action = "store",
                default=NA,
                help="number of nearest neighbour samples for SNN graph",
                metavar="numeric"),

    make_option(c("--prune"),
                type="numeric",
                action = "store",
                default=NA,
                help="prune cutoff for sample SNN graph",
                metavar="numeric"),

    make_option(c("--resolution"),
                type="numeric",
                action = "store",
                default=NA,
                help="Resolutions for louvain and leiden algorithm, numbers should be separated by comma",
                metavar="numeric"),

    make_option(c("--usegap"),
                type="logical",
                action = "store",
                default=NA,
                help="Whether to use eigengap or not",
                metavar="logical"),

    make_option(c("--sim"),
                type="logical",
                action = "store",
                default=FALSE,
                help="Is the dataset a simulated one or not",
                metavar="logical"),

    make_option(c("--truth"),
                type="character",
                action = "store",
                default='truth',
                help="Name of Column which defines ground truth of sample clusters in colData(sce)",
                metavar="character"),

    make_option(c("--nclust"),
                type="numeric",
                action = "store",
                default=NULL,
                help="Assigning number of clusters for kmeans/skmeans",
                metavar="numeric"),

    make_option(c("--graph_select"),
                type="logical",
                action = "store",
                default=NA,
                help="Whether genes should be selected on the graph",
                metavar="logical"),
    
	make_option(c("--graph_select_by_prop"),
                type="logical",
                action = "store",
                default=FALSE,
                help="Whether top variable genes should be selected by 80% criterion",
                metavar="logical"),

    make_option(c("--gcKNN"),
                type="logical",
                action = "store",
                default=NA,
                help="Whether gcKNN should be calculated",
                metavar="logical"),

    make_option(c("--SNN_mode"),
                type="character",
                action = "store",
                default=NA,
                help="SNN mode for caclust",
                metavar="character"),

    make_option(c("--overlap"),
                type="numeric",
                action = "store",
                default=NA,
                help="Overlap for graph pruning",
                metavar="numeric")
    
    

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop("Argument --file is missing.", call.=FALSE)
} else if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("Argument --outdir is missing.", call.=FALSE)
} else if (is.null(opt$name)){
    print_help(opt_parser)
    stop("Argument --name is missing.", call.=FALSE)
}


# Misc
filepath = opt$file
ntop = opt$ntop
distance = opt$distance
outdir = opt$outdir
name = opt$name
nclust = opt$nclust
sim = opt$sim
vst = opt$vst
truth = opt$truth

DECOMP <- opt$decomp


# if(nclust == "NULL"){
#     nclust <- NULL
# }

dataset = opt$dataset


# CAclust

dims = opt$dims
NNs = opt$NNs
prune = opt$prune
resol = as.numeric(opt$resolution)
usegap = opt$usegap
SNN_mode = opt$SNN_mode
graph_select <- opt$graph_select
graph_select_by_prop <- opt$graph_select_by_prop
gcKNN <- opt$gcKNN
overlap <- opt$overlap

if (is.character(overlap)) overlap <- NA


if (isTRUE(sim)){
    sim_params <- stringr::str_match(string = opt$file,
    # pattern = ".*dePROB-(?<dePROB>[0-9]_[0-9]*)_defacLOC-(?<defacLOC>[0-9]_?[0-9]*)_defacSCALE-(?<defacSCALE>[0-9]_?[0-9]*).rds")
    pattern = ".*dePROB-(?<dePROB>[0-9]_[0-9]*)_defacLOC-(?<defacLOC>[0-9]_?[0-9]*)_defacSCALE-(?<defacSCALE>[0-9]_?[0-9]*)[:graph:]{0,9}.rds")

    opt$dePROB <- as.numeric(gsub("_", ".", sim_params[,"dePROB"]))
    opt$defacLOC <- as.numeric(gsub("_", ".", sim_params[,"defacLOC"]))
    opt$defacSCALE <- as.numeric(gsub("_", ".", sim_params[,"defacSCALE"]))
}

fileformat = tools::file_ext(filepath)


if (fileformat == 'txt'){

    stop("Provided txt file as input. RDS required.")

    # cnts = read.table(opt$file, row.names = 1, header=T, sep = '\t')
    # cnts = as.matrix(cnts)
    # data = cnts

}else if (fileformat %in% c('rds', 'RDS') ){


    data = readRDS(filepath)

    if (is(data, "Seurat")){
        stop("Please provide SingleCellExperiment data, not Seurat.")
    }

    if (!is.na(ntop)){

        genevars <- modelGeneVar(data, assay.type = "logcounts")

        if (isTRUE(graph_select_by_prop) & isTRUE(graph_select)){

            chosen <- getTopHVGs(genevars, prop = 0.8, var.threshold = NULL)

        } else {

            chosen <- getTopHVGs(genevars, n = ntop, var.threshold = NULL)

        }

        data_old <- data
        data <- data[chosen,]
    }

    cnts <- as.matrix(logcounts(data))

    # genes_detect <- rowSums(cnts > 0) > (ncol(cnts)*0.01)
    # cnts <- cnts[genes_detect,]
    # data <- data[genes_detect,]

    trueclusters = colData(data)[,colnames(colData(data)) == truth]

    # if (is.null(nclust)){
    #     nclust = length(unique(colData(data)[,colnames(colData(data)) == truth]))
    # }

}


ngene = nrow(cnts)
ncell = ncol(cnts)

## write the cell indexes of each cluster to a _gt.txt file
clustnms = unique(trueclusters)
gt_file = file(file.path(outdir, paste0(dataset, '_Ntop_', ntop, '_gt.txt' )), open = 'w')

for (i in clustnms){
    idx = which(data[[truth]] == i) -1 ## substracte by 1 for cpp counting
    writeLines(paste(idx, collapse = ' '), gt_file)
}

close(gt_file)

matfile = file(file.path(outdir, paste0(dataset, '_Ntop_', ntop, '_matrix.txt' )), open = 'w')
writeLines(paste(dim(cnts), collapse = ' '), matfile)

for (i in 1:nrow(cnts)){
    writeLines(paste(rownames(cnts)[i], paste(cnts[i,], collapse = ','), sep = ',' ), matfile)
}
close(matfile)


cat("\nFinished benchmarking!\n")

