#!/bin/bash


outdir="../../discussed_data/raw/"

wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189981/suppl/GSE189981%5Fd50%5Forganoids%5Fcnt%5Fmatrix%2Etxt%2Egz" -P $outdir
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189981/suppl/GSE189981%5Fd50%5Forganoids%5Fmetadata%2Etxt%2Egz" -P $outdir

python "./convert_to_h5ad.py"

Rscript-4.2.1 "./convert_h5ad_to_sce.R"

