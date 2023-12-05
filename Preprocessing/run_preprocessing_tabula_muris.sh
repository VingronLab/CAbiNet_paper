#!/bin/bash


SCRIPT="./data_preprocessing.R"
outdir="../Data/discussed_data/preprocessed/"
file="../Data/discussed_data/raw/tabula_muris_RAW.rds"

Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name "tabula_muris_preprocessed" \
       --pct 0.01 \
       --truth "cell_ontology_class" \
       --org mm \
       --mt

