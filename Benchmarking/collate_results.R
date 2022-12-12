
library(tidyverse)
library(optparse)

option_list = list(
    make_option(c("--name"),
                type="character",
                action = "store",
                default="sample",
                help="Name of sample",
                metavar="character"),
    make_option(c("--indir"),
                type="character",
                action = "store",
                default="sample",
                help="input directory",
                metavar="character"),
    make_option(c("--outdir"),
                type="character",
                action = "store",
                default="sample",
                help="Output directory",
                metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$name)){
    print_help(opt_parser)
    stop("Argument --name is missing.", call.=FALSE)
} else if (is.null(opt$indir)){
    print_help(opt_parser)
    stop("Argument --indir is missing.", call.=FALSE)
} else if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("Argument --outdir is missing.", call.=FALSE)
}


indir <- opt$indir
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE)

runs <- list.files(path = indir,
                   pattern = "*EVALUATION.csv",
                  full.names = TRUE)

df_list <- vector(mode = "list", length = length(runs))

for (r in seq_along(runs)){

  df_list[[r]] <- read_csv(runs[r], show_col_types = FALSE, progress = FALSE)
  if ("nclust" %in% colnames(df_list[[r]])){
    if (df_list[[r]]$nclust == "NULL" | is.null(df_list[[r]]$nclust)){
        df_list[[r]]$nclust <- NA
    }
  }


}

# fills cols with NA if not present.
df_eval <- bind_rows(df_list)


saveRDS(df_eval, paste0(outdir, "/", opt$name, ".rds"))

###############################
# Example call for the script #
###############################

# Rscript collate_results.R --name "pbmc3k_benchmarking" \
#                           --indir "./results/simulated/out/pbmc3k/" \
#                           --outdir "./results/simulated/collated_results/"

