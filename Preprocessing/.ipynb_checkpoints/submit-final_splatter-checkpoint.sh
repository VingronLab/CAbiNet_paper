
THREADS=4
MEMORY=50G
MINUTES=120
TMPDIR=10G

indir="sim_data/"

outdir="sim_data/preprocessed/"
logdir=${outdir}_log
SCRIPT="data_preprocessing_splatter.R"

mkdir -p $outdir
# cd $outdir

mkdir -p $logdir

# files=`readlink -f /project/CAclust/data/CAbiNET_benchmarking_data/data/*.rds`
truth='Group'
pcts=0.01
ntop=NULL


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