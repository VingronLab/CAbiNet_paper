#!/bin/bash

THREADS=16
MEMORY=8G
MINUTES=600

indir=./results/out
outdir=./results/collate_res
mkdir -p $outdir

name='sim_2kgenes'

 mxqsub --stdout="${outdir}/$name.stdout.log" \
       --group-name=$name \
       --threads=$THREADS \
       --memory=$MEMORY \
       -t $MINUTES \
       Rscript-4.2.1 collate_results.R --indir $indir \
                                       --outdir $outdir \
                                       --name $name
