#!/bin/bash


#########################################################################################################
## Parameters for mxq-cluster job submission, change it according to your in-house server requirements.##
#########################################################################################################

THREADS=4
MEMORY=100G
MINUTES=120
TMPDIR=80G


# indir="./ExperimentalData/rawdata"

indir="../Data/rawdata/realdata"
outdir="./ExperimentalData/preprocessed"
logdir=${outdir}_log
SCRIPT="./data_preprocessing.R"

mkdir -p $outdir
# cd $outdir

mkdir -p $logdir

truth='truth'
pcts=0.01
ntop=NULL


###############
# ZeiselBrain #
###############

file="${indir}/ZeiselBrain.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org mm \
       --mt

#################
# BaronPancreas #
#################

file="${indir}/BaronPancreas.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt


############
# Darmanis #
############

file="${indir}/Darmanis.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt

###############
# FreytagGold #
###############

file="${indir}/FreytagGold.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt

################
# HeOrganAtlas #
################

file="${indir}/HeOrganAtlas.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs

############
# PBMC_10X #
############

file="${indir}/PBMC_10X.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt

########################
# Tirosh_nonmaglignant #
########################

file="${indir}/Tirosh_nonmaglignant.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt

###################
# brain_organoids #
###################

file="${indir}/brain_organoids.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs


#########################
# tabula_sapiens_tissue #
#########################

file="${indir}/tabula_sapiens_tissue.rds"
filename=`basename $file .rds`

mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
       --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
       --group-name="${filename}" \
       --threads=$THREADS \
       --memory=$MEMORY \
       --tmpdir=$TMPDIR \
       -t $MINUTES \
      Rscript-4.2.1 $SCRIPT   \
       --outdir $outdir \
       --file $file \
       --name $filename \
       --pct $pcts \
       --truth $truth \
       --org hs \
       --mt

 #########################
 # dnel spatial E14-16h #
 #########################

 file="${indir}/dmel_E14-16h.rds"
 filename=`basename $file .rds`

 mxqsub --stdout="${logdir}/${filename}_filtered${pcts}.stdout.log" \
        --stderr="${logdir}/${filename}_filtered${pcts}.stderr.log" \
        --group-name="${filename}" \
        --threads=$THREADS \
        --memory=$MEMORY \
        --tmpdir=$TMPDIR \
        -t $MINUTES \
       Rscript-4.2.1 $SCRIPT   \
        --outdir $outdir \
        --file $file \
        --name $filename \
        --pct $pcts \
        --truth $truth \
        --org hs \
        --rmbatch TRUE \
        --batchlab 'slice_ID' \
        --modpoisson FALSE \
        --mt


# for f in $files; do

#     for p in ${pcts[@]}; do

#         filename=`basename $f .rds`

#         if [ $f = "ZeiselBrain.rds" ]; then

#             mxqsub --stdout="${logdir}/${filename}_${m}.stdout.log" \
#                    --stderr="${logdir}/${filename}_${ncell}.stderr.log" \
#                    --group-name="${filename}" \
#                    --threads=$THREADS \
#                    --memory=$MEMORY \
#                    --tmpdir=80G \
#                    -t $MINUTES \
#                   Rscript-4.2.1 $SCRIPT   \
# 		           --outdir $outdir \
#                    --file $f \
#                    --name $filename \
#                    --pct $p \
#                    --org mm
#         else
#             mxqsub --stdout="${logdir}/${filename}/${filename}_${m}.stdout.log" \
#                    --stderr="${logdir}/${filename}_${ncell}_pct${p}.stderr.log" \
#                    --group-name="${filename}" \
#                    --threads=$THREADS \
#                    --memory=$MEMORY \
#                    --tmpdir=80G \
#                    -t $MINUTES \
#                   Rscript-4.2.1 $SCRIPT   \
# 		           --outdir $outdir \
#                    --file $f \
#                    --name $filename \
#                    --pct $p \
#                    --truth $truth \
#                    --org hs

#         fi
#     done
# done
