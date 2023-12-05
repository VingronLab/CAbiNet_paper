
library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(scDblFinder)
library(Seurat)
library(optparse)

set.seed(2358)


option_list = list(
    make_option(c("--name"),
                type="character",
                action = "store",
                default = NULL,
                help="Name of sample",
                metavar="character"),

    make_option(c("--file"),
                type="character",
                action = "store",
                default = NULL,
                help="Path to file",
                metavar="character"),

    make_option(c("--outdir"),
                type="character",
                action = "store",
                default = NULL,
                help="Output directory",
                metavar="character"),

    make_option(c("--org"),
                type="character",
                action = "store",
                default = NULL,
                help="Organism. hs or mm",
                metavar="character"),
    make_option(c("--mt"),
                type="logical",
                action = "store_true",
                default=FALSE,
                help="Whether to use mt filtering or not",
                metavar="logical"),
      make_option(c("--truth"),
                type="character",
                action = "store",
                default = "phenoid",
                help="the name of column in metadata which gives ground-truth of celltypes",
                metavar="character"),
      make_option(c("--rmbatch"),
                type="logical",
                action = "store",
                default = FALSE,
                help="Logical, remove batch effect or not",
                metavar="logiacal"),
      make_option(c("--batchlab"),
                type="character",
                action = "store",
                default = NULL,
                help="name of column where the batch information is stored in colData(sce)",
                metavar="character"),
      make_option(c("--modpoisson"),
                type="logical",
                action = "store",
                default = TRUE,
                help="whether model variance of genes by Possion dist or not",
                metavar="logical"),

    make_option(c("--pct"),
                type="numeric",
                action = "store",
                default = NULL,
                help="genes expressed in less than 'pct'% will be filtered out",
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
} else if (is.null(opt$org)){
    print_help(opt_parser)
    stop("Argument --org is missing.", call.=FALSE)
}

name <- opt$name
file <- opt$file
outdir <- opt$outdir
org <- opt$org
mt_filter <- opt$mt
pct <- opt$pct
truth <- opt$truth
rmbatch <- opt$rmbatch
batchlab <- opt$batchlab
modPoisson <- opt$modpoisson

if (is.null(name)) name <- tools::file_path_sans_ext(basename(file))

outdir <- file.path(outdir, name)
imgdir <- file.path(outdir, "img")

dir.create(outdir)
dir.create(imgdir)

#######################

sce <- readRDS(file)

if(is(sce, 'Seurat')){
    sce <- SingleCellExperiment(list(counts = GetAssayData(sce,
                                                           slot = 'count',
                                                           assay = 'RNA')),
                                colData = sce@meta.data)
}


# Remove genes not expressed in any cells
sce <- sce[rowSums(counts(sce)) > 0,]


###########################################3

# change gene names to gene symbol
getSCEWithSymbols <- function(sce, keytype = "ENTREZID", org.db = org.Mm.eg.db){
    if (keytype == 'SYMBOL'){
        rowData(sce)$SYMBOL <- rownames(sce)
        return(sce)
    }
    require("AnnotationDbi")
    require("org.Hs.eg.db")
    require("org.Mm.eg.db")

    geneSymbols <- mapIds(org.db, keys=rownames(sce), column="SYMBOL", keytype= keytype, multiVals="first")

    rowData(sce)$SYMBOL <- geneSymbols
    return(sce)
}

template = rownames(sce)[1]

if (!isEmpty(grep('EN', template))){

    names <- stringr::word(rownames(sce), 1, 1, '_')
    rownames(sce) <- names

    if (!isEmpty(grep('ENSG', template))){
        org.db <- org.Hs.eg.db
    }

    if (!isEmpty(grep('ENSMU', template))){
        org.db <- org.Mm.eg.db
        org <- 'mm'
    }

    if (!isEmpty(grep('[.]', template))){

        rownames(sce) <- stringr::word(rownames(sce), 1, 1, '[.]')

    }

    sce <- getSCEWithSymbols(sce, keytype = 'ENSEMBL', org.db = org.db)
} else {
    rowData(sce)$SYMBOL <- rownames(sce)
}


# Filter based on QC metrics

if(org == "hs"){
    mt_genes <- grepl("^MT-",  rowData(sce)$SYMBOL)

} else if (org == "mm"){

  mt_genes <- grepl("^mt-",  rowData(sce)$SYMBOL)

}

if(isTRUE(mt_filter)){

    qc_df <- perCellQCMetrics(sce, subsets=list(Mito=mt_genes))

    reasons <- perCellQCFilters(qc_df,
                                sum.field = "sum",
                                detected.field = "detected",
                                sub.fields=c("subsets_Mito_percent"))
} else {

    qc_df <- perCellQCMetrics(sce, subsets=list(Mito=mt_genes))

    reasons <- perCellQCFilters(qc_df,
                                sum.field = "sum",
                                detected.field = "detected",
                                sub.fields = NULL)
}



# QC plots

colData(sce) <- cbind(colData(sce), qc_df)
sce$discard <- reasons$discard

p_counts <- plotColData(sce, x=truth, y="sum", colour_by="discard") +
                scale_y_log10() +
                ggtitle("Total count") +
                theme(axis.text.x = element_text(angle = 90))

ggsave(p_counts, file = file.path(imgdir, paste(name, "total_counts_violinplot.png", sep = "_")))

p_detected <- plotColData(sce, x=truth, y="detected", colour_by="discard") +
              scale_y_log10() +
              ggtitle("Detected features") +
              theme(axis.text.x = element_text(angle = 90))

ggsave(p_detected, file = file.path(imgdir, paste(name, "detected_violinplot.png", sep = "_")))

p_mito <- plotColData(sce, x=truth, y="subsets_Mito_percent", colour_by="discard") +
              ggtitle("Mito percent") +
              theme(axis.text.x = element_text(angle = 90))

ggsave(p_mito, file = file.path(imgdir, paste(name, "percent_mito_violinplot.png", sep = "_")))


p_sumVSdetected <- plotColData(sce, x="sum", y="detected", colour_by="discard") +
              ggtitle("sum vs detected")

ggsave(p_sumVSdetected, file = file.path(imgdir, paste(name, "sum_vs_detected_scatter.png", sep = "_")))

p_hist_detected <- ggplot(as.data.frame(colData(sce)), aes(x=detected))+
  geom_histogram(fill = "#008080") +
  ggtitle("detected") +
  theme_bw()

ggsave(p_hist_detected, file = file.path(imgdir, paste(name, "detected_hist.png", sep = "_")))

p_hist_sum <- ggplot(as.data.frame(colData(sce)), aes(x=sum))+
  geom_histogram(fill = "#008080") +
  ggtitle("sum") +
  theme_bw()

ggsave(p_hist_sum, file = file.path(imgdir, paste(name, "sum_hist.png", sep = "_")))



# remove bad cells
sce <- sce[,!reasons$discard]

# filter out genes expressed in less than 1% cells
bm <- counts(sce) > 0

ncell <- round(ncol(sce) * pct)
label <- paste0('pct', pct)

idx <- which((rowSums(bm) > ncell))
sce <- sce[idx, ]


# Normalization
clust.sce <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)
sce <- logNormCounts(sce)

if (isTRUE(rmbatch)){
  if (is.null(batchlab)){
    stop('value of parameter batchlab is missing.')
  }
  cnt_bc = sva::ComBat(counts(sce), batch = colData(sce)[, batchlab])
  if (min(cnt_bc) < 0){
    cnt_bc <- cnt_bc - min(cnt_bc)
  }
  cnt_bc <- cnt_bc[rowSums(cnt_bc) > 0 , colSums(cnt_bc) > 0]
  counts(sce) <- cnt_bc
  sce <- logNormCounts(sce)
}


#dim reduc.
if (isTRUE(modPoisson)){
  sce.dec <- modelGeneVarByPoisson(sce)
  sce.top <- getTopHVGs(sce.dec, prop = 0.2)
}else{
  sce.dec <- modelGeneVar(sce)
  sce.top <- getTopHVGs(sce.dec, prop = 0.2)
}


sce <- denoisePCA(sce, technical=sce.dec, subset.row=sce.top)
sce <- runUMAP(sce, dimred="PCA")

p_umap<- plotUMAP(sce, colour_by=truth)
ggsave(p_umap, file = file.path(imgdir, paste(name, "Truth_umap.png", sep = "_")))

# doublet detection

sce <- scDblFinder(sce, clusters=clust.sce)
# table(sce$scDblFinder.class)

p_umap_dblscore <- plotUMAP(sce, colour_by="scDblFinder.score")
ggsave(p_umap_dblscore, file = file.path(imgdir, paste(name, "doublet_score_umap.png", sep = "_")))

p_umap_dblclass <- plotUMAP(sce, colour_by="scDblFinder.class")
ggsave(p_umap_dblclass, file = file.path(imgdir, paste(name, "doublet_class_umap.png", sep = "_")))

# sce <- sce[,sce$scDblFinder.class != "doublet"]




# save preprocessed file
saveRDS(sce, file.path(outdir, paste0(name, "_filtered.rds")))
