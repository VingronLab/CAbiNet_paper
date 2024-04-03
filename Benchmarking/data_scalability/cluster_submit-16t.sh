#!/bin/bash

#########################################################################################################
########### Attention: 
########### The simulated data sets can be downloaded from online data repository.
########### Or you can also run the sim_generation.R script to generate simulated data sets by yourself.
#########################################################################################################

# Monocle3 requires specific libraries.
LD_LIBRARY_VAR=/project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/renv/local/lib:../renv/local/include
# LD_LIBRARY_VAR=$HOME/.local/lib:$HOME/.local/bin:$HOME/.local/include
GDAL_DATA_VAR=/project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/renv/local/share/gdal
# GDAL_DATA_VAR=$HOME/.local/share/gdal


THREADS=16
MEMORY=1000G
MINUTES=600

cd ../
#scripts_path="../algorithms"
scripts_path='data_scalability/algorithms-parallel'

outdir=./data_scalability/results_16t_parallel
mkdir -p $outdir

logdir="${outdir}/log/"
mkdir -p $outdir
mkdir -p $logdir

OUTDIR="${outdir}/out/"
mkdir -p $OUTDIR

here_dir="${outdir}/sh/"
mkdir -p $here_dir
mkdir -p $here_dir/.done/

#indir='../../Data/sim_data/raw/data_scalability'
indir='/project/CAclust/CAbiNET_paper/paper_revision/benchmarking/data_scalability/pbmc3k'

datasets=`ls -d ${indir}/*`

for dataset in ${datasets[@]}; do

dataset=`basename $dataset`

    
files="${indir}/${dataset}/*.rds"

ntop=(2000)
nclust=6
truth='Group'


####################################################################################
# CAclust 
res=1
dims=(20)
NNs=(20)
overlap=(NA)
graph_select=(TRUE)
SNN_mode=("out")
calc_gckNN=(FALSE)
prune=$(echo 'print(1/15)' | python3)





for f in ${files[@]}; do

    filename=`basename $f .rds`

    for nt in ${ntop[@]}; do



    ###########
    # CAclust #
    ###########

    for d in ${dims[@]}; do

        for n in ${NNs[@]}; do

            for gs in ${graph_select[@]}; do

                for snn in ${SNN_mode[@]}; do

                    for gc in ${calc_gckNN[@]}; do

                        counter=0

                        for ov in ${overlap[@]}; do

                        if [ $gs = "FALSE" ]; then

                            if (($counter >= 1)); then
                                break
                            else
                                ov=NA
                                MEM_caclust=500G
                            fi

                        else
                            MEM_caclust=$MEMORY
                        fi
                        
                        echo $MEM_caclust

                        algorithm="caclust_leiden"
                        SCRIPT="${scripts_path}/${algorithm}.R"

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"

tmp_sh="${here_dir}/bench_sim_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU
cd /project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/Benchmarking

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --dims $d    \\
   --prune $prune \\
   --NNs $n     \\
   --resolution $res     \\
   --sim TRUE \\
   --truth $truth \\
   --graph_select $gs \\
   --prune_overlap FALSE \\
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                        mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
                               --group-name="bench_sim_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh
 
                       algorithm="caclust_igraph"
                        SCRIPT="${scripts_path}/${algorithm}.R"
echo $algorithm

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"

tmp_sh="${here_dir}/bench_sim_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU
cd /project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/Benchmarking

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --dims $d    \\
   --prune $prune \\
   --NNs $n     \\
   --resolution $res     \\
   --sim TRUE \\
   --truth $truth \\
   --graph_select $gs \\
   --prune_overlap FALSE \\
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                        mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
                               --group-name="bench_sim_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh
<<!
                        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_sim_scalability_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               Rscript-4.2.1 $SCRIPT   \
                                   --outdir $OUTDIR  \
                                   --file $f \
                                   --dataset $dataset \
                                   --name $nm \
                                   --ntop $nt \
                                   --dims $d    \
                                   --prune $prune \
                                   --prune_overlap FALSE \
                                   --NNs $n     \
                                   --resolution $res     \
                                   --sim TRUE \
                                   --truth $truth \
                                   --graph_select $gs \
                                   --nclust NULL \
                                   --SNN_mode $snn \
                                   --gcKNN $gc \
                                   --overlap $ov
!
                        algorithm="caclust_spectral"
                        SCRIPT="${scripts_path}/${algorithm}.R"

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"

tmp_sh="${here_dir}/bench_sim_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU
cd /project/CAclust_scripts/CAclust_paper/2024_NAR_submission_files/Benchmarking

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --dims $d    \\
   --prune $prune \\
   --NNs $n     \\
   --resolution $res     \\
   --sim TRUE \\
   --truth $truth \\
   --graph_select $gs \\
   --usegap TRUE \\
   --prune_overlap FALSE \\
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                        mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
                               --group-name="bench_sim_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh

<<!
                        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_sim_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               Rscript-4.2.1 $SCRIPT   \
                                   --outdir $OUTDIR  \
                                   --file $f \
                                   --dataset $dataset \
                                   --name $nm \
                                   --ntop $nt \
                                   --dims $d    \
                                   --prune $prune \
                                   --prune_overlap FALSE \
                                   --NNs $n     \
                                   --resolution $res     \
                                   --sim TRUE \
                                   --truth $truth \
                                   --usegap TRUE \
                                   --graph_select $gs \
                                   --nclust NULL \
                                   --SNN_mode $snn \
                                   --gcKNN $gc \
                                   --overlap $ov
!
                        counter=$((counter+1))
                        
                        done
                    done
                done                              
            done
        done
    done


    # end caclust


    done
done
done
