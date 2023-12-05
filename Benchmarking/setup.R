library(Rcpp)
library(CAbiNet)
library(APL)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
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
library(dplyr)
library(monocle3)


set.seed(2358)

source("./sim_eval.R")
source("./algorithms/biclustlib/clustering_error.R")


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

    make_option(c("--prune_overlap"),
                type="logical",
                action = "store",
                default=TRUE,
                help="prune gene nodes or not in the SNN graph by overlapping of neighbourhood",
                metavar="logical"),

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
    
    make_option(c("--leiden_pack"),
                type="character",
                action = "store",
                default='igraph',
                help="package for running leiden algorithm",
                metavar="character"),

    make_option(c("--overlap"),
                type="numeric",
                action = "store",
                default=NA,
                help="Overlap for graph pruning",
                metavar="numeric"),


    # QUBIC options
    make_option(c("--r_param"),
                type="numeric",
                action = "store",
                default=NA, # default 1
                help="The range of possible ranks",
                metavar="numeric"),

    make_option(c("--q_param"),
            type="numeric",
            action = "store",
            default=NA, # default 0.06
            help="QUBIC q parameter",
            metavar="numeric"),

    make_option(c("--c_param"),
            type="numeric",
            action = "store",
            default=NA, # default 0.95
            help="QUBIC c param",
            metavar="numeric"),
    #s4vd parameters

    make_option(c("--pcerv"),
            type="numeric",
            action = "store",
            default=NA, # default 0.05
            help="Per comparsion wise error rate for v.",
            metavar="numeric"),

    make_option(c("--pceru"),
        type="numeric",
        action = "store",
        default=NA, # default 0.05
        help="Per comparsion wise error rate for u.",
        metavar="numeric"),

    make_option(c("--ss_thr_min"),
        type="numeric",
        action = "store",
        default=NA, # default 0.6
        help="Range of the cutoff threshold minimum.",
        metavar="numeric"),

    make_option(c("--ss_thr_add"),
        type="numeric",
        action = "store",
        default= NA, # default 0.05
        help="Range of the cutoff threshold: to add on minimum",
        metavar="numeric"),
    #Plaid parameters

    make_option(c("--rrelease"),
        type="numeric",
        action = "store",
        default= NA, # default 0.7
        help="threshold to prune rows in the layers depending on row homogeneity",
        metavar="numeric"),

    make_option(c("--crelease"),
        type="numeric",
        action = "store",
        default= NA, # default 0.7
        help="threshold to prune rows in the layers depending on column homogeneity",
        metavar="numeric"),
    #Unibic parameters

    make_option(c("--t_param"),
        type="numeric",
        action = "store",
        default= NA, # default 0.95
        help="consistency level of the block (0.5-1.0]",
        metavar="numeric"),

    make_option(c("--q_discr"),
        type="numeric",
        action = "store",
        default= NA, # default 0
        help="a double value for quantile discretization",
        metavar="numeric"),

    # Bimax parameters
    make_option(c("--minr"),
        type="numeric",
        action = "store",
        default= NA, # default 2
        help="Minimum row size of resulting bicluster",
        metavar="numeric"),

    make_option(c("--minc"),
        type="numeric",
        action = "store",
        default= NA, # default 2
        help="Minimum column size of resulting bicluster",
        metavar="numeric"),
    # make_option(c("--maxc"),
    #     type="numeric",
    #     action = "store",
    #     default= NA, # default 12
    #     help="Maximum column size of resulting bicluster",
    #     metavar="numeric"),
    # CCA

    make_option(c("--alpha"),
        type="numeric",
        action = "store",
        default= NA, # default 1.5
        help="Scaling factor",
        metavar="numeric"),

    make_option(c("--delta"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="Maximum of accepted score",
        metavar="numeric"),

    # Xmotifs
    make_option(c("--ns_param"),
        type="numeric",
        action = "store",
        default= NA, # default 10
        help="Number of columns choosen",
        metavar="numeric"),

    make_option(c("--nd_param"),
        type="numeric",
        action = "store",
        default= NA, # default 10
        help="Number of repetitions.",
        metavar="numeric"),

    make_option(c("--sd_param"),
        type="numeric",
        action = "store",
        default= NA, # default 5
        help="Sample size in repetitions.",
        metavar="numeric"),

    make_option(c("--alphaX"),
        type="numeric",
        action = "store",
        default= NA, # default 0.05
        help="Scaling factor for column result.",
        metavar="numeric"),
    # BC_Spectral

    make_option(c("--numEig"),
        type="numeric",
        action = "store",
        default= NA, # default 1.5
        help="the number of eigenValues considered to find biclusters.",
        metavar="numeric"),

    make_option(c("--minr_BCS"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="minimum number of rows that biclusters must have.",
        metavar="numeric"),

    make_option(c("--minc_BCS"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="Sminimum number of columns that biclusters must have.",
        metavar="numeric"),

    make_option(c("--withinVar"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="maximum within variation allowed.",
        metavar="numeric"),

    # make_option(c("--nbest"),
    #     type="numeric",
    #     action = "store",
    #     default= NA, # default 1
    #     help="number of eigenvectors to which the data is projected for the final clustering step",
    #     metavar="numeric"),

    # make_option(c("--nclusters"),
    #     type="numeric",
    #     action = "store",
    #     default= NA, # default 1
    #     help="vector with first element the number of row clusters and second element the number of column clusters.",
    #     metavar="numeric")

    make_option(c("--Qconsistency"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="consistency level of the block (0.5-1.0],",
        metavar="numeric"),

    make_option(c("--Qoverlap"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="filtering overlapping blocks",
        metavar="numeric"),

    make_option(c("--Qcmin"),
        type="numeric",
        action = "store",
        default= NA, # default 1
        help="minimum column width of the block.",
        metavar="numeric"),

    # Seurat options
    make_option(c("--logfc_thr"),
        type="numeric",
        action = "store",
        default= NA, # default 0.25
        help="minimum log2fc threshold to test genes.",
        metavar="numeric"),
    make_option(c("--min_perc"),
        type="numeric",
        action = "store",
        default= NA, # default 0.25
        help="minimum fraction of min.pct cells to test genes in.",
        metavar="numeric"),
    make_option(c("--rthr"),
        type="numeric",
        action = "store",
        default= NA, # default 0.01
        help="return threshold above which genes are returned.",
        metavar="numeric"),

    # Monocle3
    make_option(c("--redm"),
        type="character",
        action = "store",
        default= 'UMAP',
        help="Dimensin reduction method of Monocle3",
        metavar="character"),
    make_option(c("--ngene_pg"),
        type="numeric",
        action = "store",
        default= 100, # default 0.25
        help="Number of marker genes to calculate by Monocle3",
        metavar="numeric"),
    
    # DivBiclust
    make_option(c("--maxdiff"),
        type="numeric",
        action = "store",
        default= 0.15,
        help="Argument max_diff",
        metavar="numeric"),
    make_option(c("--dorate"),
        type="numeric",
        action = "store",
        default= 0.1,
        help="fraction of missing values in a bicluster",
        metavar="numeric"),
    make_option(c("--seedColSz"),
        type="numeric",
        action = "store",
        default= 50,
        help="size of seed gene set",
        metavar="numeric"),
     make_option(c("--maxColSz"),
        type="numeric",
        action = "store",
        default= 100,
        help="maximum size of gene set, fixed to 100",
        metavar="numeric"),
    make_option(c("--simThresh"),
        type="numeric",
        action = "store",
        default= 0.5,
        help="similarity threshold for pattern merging, fixed to 0.5",
        metavar="numeric"),
    make_option(c("--isdivbiclust"),
         type = 'logical',
         action = 'store',
         default = FALSE,
         help = 'whether the chosen algorithm is divbiclust or not',
         metavar = 'logical'),

    # BackSpin
    
    make_option(c("--numLevels"),
        type="numeric",
        action = "store",
        default = 2,
        help="the number of splits that will be tried",
        metavar="numeric"),
    
    make_option(c("--stop_const"),
        type="numeric",
        action = "store",
        default= 1.15,
        help="minimum score that a breaking point has to reach to be suitable for splitting",
        metavar="numeric"),
    
    make_option(c("--low_thrs"),
        type="numeric",
        action = "store",
        default= 0.2,
        help="genes with average lower than this threshold are assigned to either of the splitting group reling on genes that are higly correlated with them",
        metavar="numeric")
    
    #make_option(c("--pcapreproc"),
    #     type = 'logical',
    #     action = 'store',
    #     default = FALSE,
    #     help = 'perform PCA dim. reduction as a preprocessing step.',
    #     metavar = 'logical')

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

dataset=opt$dataset


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

#pca_preproc <- opt$pcapreproc

if(nclust == "NULL"){
    nclust <- NULL
}

# CAclust

dims = opt$dims
NNs = opt$NNs
prune = opt$prune
prune_overlap = opt$prune_overlap
resol = as.numeric(opt$resolution)
usegap = opt$usegap
SNN_mode = opt$SNN_mode
graph_select <- opt$graph_select
graph_select_by_prop <- opt$graph_select_by_prop
gcKNN <- opt$gcKNN
overlap <- opt$overlap
leiden_pack = opt$leiden_pack

if (is.character(overlap)) overlap <- NA

# QUBIC
c_param <- opt$c_param
r_param <- opt$r_param
q_param <- opt$q_param

#s4vd

pcerv <- opt$pcerv
pceru <- opt$pceru
ss_thr <- c(opt$ss_thr_min, (opt$ss_thr_min+opt$ss_thr_add))

# Plaid parameters
rrelease <- opt$rrelease
crelease <- opt$crelease

# Unibic
t_param <- opt$t_param
q_discr <- opt$q_discr

# Bimax
minr <- opt$minr
minc <- opt$minc
# maxc <- opt$maxc

# CCA
alpha <- opt$alpha
delta <- opt$delta

# Xmotifs
ns_param <- opt$ns_param
nd_param <- opt$nd_param
sd_param <- opt$sd_param
alphaX <- opt$alphaX

# BC_Spectral
numEig <- opt$numEig
minr_BCS <- opt$minr_BCS
minc_BCS <- opt$minc_BCS
withinVar <- opt$withinVar
nbest <- opt$nbest
nclusters <- opt$nclusters

if (!is.numeric(nclusters)){
    nclusters <- NULL
}

# IRISFGM/QUBIC2

Qconsistency <- opt$Qconsistency
Qoverlap <- opt$Qoverlap
Qcmin <- opt$Qcmin

# Seurat
logfc_thr <- opt$logfc_thr
min_perc <- opt$min_perc
rthr <- opt$rthr

# Monocle
ntop = opt$ntop
resolution = as.numeric(opt$resolution)
reduction_method = opt$redm
genes_to_test_per_group = opt$ngene_pg

# DivBiclust
max_diff = opt$maxdiff
do_rate = opt$dorate
seed_col_sz = opt$seedColSz
max_col_sz = opt$maxColSz
simThresh = opt$simThresh
is_divbiclust = opt$isdivbiclust


# BackSpin
numL <- opt$numLevels
stopc <- opt$stop_const
lowT <- opt$low_thrs


if (isTRUE(sim)){
    sim_params <- stringr::str_match(string = opt$file,
    # pattern = ".*dePROB-(?<dePROB>[0-9]_[0-9]*)_defacLOC-(?<defacLOC>[0-9]_?[0-9]*)_defacSCALE-(?<defacSCALE>[0-9]_?[0-9]*).rds")
    pattern = ".*dePROB-(?<dePROB>[0-9]_[0-9]*)_defacLOC-(?<defacLOC>[0-9]_?[0-9]*)_defacSCALE-(?<defacSCALE>[0-9]_?[0-9]*)[:graph:]{0,9}.rds")

    opt$dePROB <- as.numeric(gsub("_", ".", sim_params[,"dePROB"]))
    opt$defacLOC <- as.numeric(gsub("_", ".", sim_params[,"defacLOC"]))
    opt$defacSCALE <- as.numeric(gsub("_", ".", sim_params[,"defacSCALE"]))
}



if (isTRUE(is_divbiclust)){
    # prefix of input file name
    ds_type = file.path(outdir,paste0(name, '_Ntop_', ntop))
   	in_file = file.path(outdir,paste0(dataset, '_Ntop_', ntop)) 
    # calculate the dropout rates in the input data set (the data that has been preprocessed)
    filepath = file.path(outdir, paste0(dataset, '_Ntop_', ntop, '_matrix.txt' ))
    mat = read.csv(filepath,  header = F,  skip = 1)
    mat = mat[, 2:ncol(mat)]
    # mat[is.na(mat)] = 0
    # do_rate = sum(mat == 0)/dim(mat)[1]/dim(mat)[2]
    ncell = ncol(mat)
                    
}

fileformat = tools::file_ext(filepath)

if (fileformat == 'txt'){

    # stop("Provided txt file as input. RDS required.")
    cat('running divbiclust.....')

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


    trueclusters = colData(data)[,colnames(colData(data)) == truth]

    if (is.null(nclust)){
        nclust = length(unique(colData(data)[,colnames(colData(data)) == truth]))
    }

}
