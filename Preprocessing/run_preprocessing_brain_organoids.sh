#!/bin/bash


SCRIPT="./data_preprocessing_brain_organoids.R"
outdir="../Data/discussed_data/preprocessed/"
file="../Data/discussed_data/raw/brain_organoids_RAW.rds"

Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name "brain_organoids_preprocessed" \
       --pct 0.01 \
       --truth "cell_type" \
       --org hs \
       --mt_perc 40 \
       --mt

