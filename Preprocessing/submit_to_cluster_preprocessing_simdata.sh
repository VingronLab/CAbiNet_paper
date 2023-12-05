#!/bin/bash


#########################################################################################################
## Parameters for mxq-cluster job submission, change it according to your in-house server requirements.##
#########################################################################################################

THREADS=4
MEMORY=50G
MINUTES=120
TMPDIR=10G

# indir="./sim_data/sim_data"

indir="../Data/sim_data/raw"
outdir="../Data/sim_data/preprocessed"

logdir=${outdir}_log
SCRIPT="./data_preprocessing_splatter.R"

mkdir -p $outdir
mkdir -p $logdir

truth='Group'
pcts=0.01 # percent of cells a gene needs to be expressed in.
ntop=NULL


###############################################
#### Sim. data based on Zeisel Brain data #####
###############################################

dataset="zeisel"

files="${indir}/${dataset}/*.rds"

resdir="${outdir}/${dataset}"
mkdir -p $resdir

for f in ${files[@]}; do

       filename=`basename $f .rds`

       mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
              --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
              --group-name="${filename}" \
              --threads=$THREADS \
              --memory=$MEMORY \
              --tmpdir=$TMPDIR \
              -t $MINUTES \
             Rscript-4.2.1 $SCRIPT   \
              --outdir $resdir \
              --file $f \
              --name $filename \
              --pct $pcts \
              --truth $truth \


done


#########################################
#### Sim. data based on pbmc3k data #####
#########################################

dataset="pbmc3k"

files="${indir}/${dataset}/*.rds"

resdir="${outdir}/${dataset}"
mkdir -p $resdir

for f in ${files[@]}; do

       filename=`basename $f .rds`

       mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
              --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
              --group-name="${filename}" \
              --threads=$THREADS \
              --memory=$MEMORY \
              --tmpdir=$TMPDIR \
              -t $MINUTES \
             Rscript-4.2.1 $SCRIPT   \
              --outdir $resdir \
              --file $f \
              --name $filename \
              --pct $pcts \
              --truth $truth \

done
