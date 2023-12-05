#!/bin/bash

THREADS=6
MEMORY=4G
MINUTES=240


#is_preprocessing=true
is_preprocessing=false


#dataset="pbmc3k"
dataset="zeisel"


scripts_path="./algorithms/"

outdir="./results/simulated/"

logdir="${outdir}/log/${dataset}"
here_dir="${outdir}/sh/${dataset}"

mkdir -p $outdir
mkdir -p $logdir
mkdir -p $here_dir
mkdir -p $here_dir/.done/

OUTDIR="${outdir}/out/${dataset}"
mkdir -p $OUTDIR

files="../Data/sim_data/preprocessed/${dataset}/*.rds"

ntop=(2000 4000 6000)
nclust=6
truth='Group'

if [ $is_preprocessing = true ] ; then

######################################################################################
######## run pre-divbiclust.r to prepare the input matrix and ground-truth labels#####
######## this will ouput two '.txt' files which will be used as input for divbiclust##
######################################################################################

OUTDIR="${outdir}/out/${dataset}"
mkdir -p $OUTDIR

for f in ${files[@]}; do

    filename=`basename $f .rds`

    for nt in ${ntop[@]}; do
   
       
            algorithm="divbiclust"
                        SCRIPT="${scripts_path}/pre-divbiclust.R"
    
                        nm="${algorithm}_${filename}_ntop-${nt}"
    
                        mxqsub --stdout="${logdir}/bench_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEMORY \
                               -t $MINUTES \
                               Rscript-4.2.1 $SCRIPT   \
                                   --outdir $OUTDIR  \
                                   --file $f \
                                   --dataset $filename \
                                   --truth $truth \
                                   --name $nm \
                                   --sim TRUE \
                                   --ntop $nt
    done
done

####################################################################################

else


# CAclust - 108
res=1
dims=(40 80 150)
NNs=(30 60 100)
overlap=(0.1)
graph_select=(FALSE TRUE)
SNN_mode=("all")
calc_gckNN=(TRUE FALSE)
prune=$(echo 'print(1/15)' | python3)

# QUBIC - 108
r_param=(1 3 6)
c_param=(0.90 0.95 0.99)
q_param=(0.01 0.06 0.1 0.2)

# s4vd - 108
pcerv=(0.01 0.05 0.15)
pceru=(0.01 0.05 0.15)
ss_thr_min=(0.6 0.7)
ss_thr_add=(0.05 0.15)


# Plaid - 108
rrelease=(0.50 0.54 0.58 0.62 0.66 0.70)
crelease=(0.50 0.54 0.58 0.62 0.66 0.70)

# Unibic - 108
t_param=(0.50 0.59 0.68 0.77 0.86 0.95)
q_discr=(0.0 0.1 0.2 0.3 0.4 0.5)


# Bimax - 108
minr=(2 40 90 140 190 240)
minc=(2 12 24 36 48 60)
# maxc=(300 600 1000)

# CCA - 108
delta=(0.5 1.0 1.5 2.0 2.5 3.0)
alpha=(0.5 1.0 1.5 2.0 2.5 3.0)

# Xmotifs - 108
ns_param=(30 150 400)
nd_param=(10)
sd_param=(5 10 20)
alpha_X=(0.01 0.05 0.1 0.5)

# BC_Spectal - 108
numEig=(3 9 15 21 27 33)
minr_BCS=(2 50 100)
minc_BCS=(2 60)
withinVar=(1)
# nbest=(2 3)
# n_clusters=(NULL ${nclust})

# IRISFGM/QUBIC2 - 108
qcons=(0.5 0.66 0.83 1)
qoverlap=(0.5 0.75 1)
qcmin=(10 50 100)


# Seurat
res_seurat=(1)
dims_seurat=(10 30 60)
NNs_seurat=(20 40 80)
min_perc=(0.15 0.25)
logfc_thr=(0.25)
return_thr=(0.01 0.05)

# Monocle 108
res_monocle=(0.1 1)
dims_monocle=(10 30 60)
NNs_monocle=(20 40 80)
reduce_method=("UMAP" "PCA")
ngene_perg=(50)

# divbiclust
maxdiffs=(0.15 0.2 0.25)
seedColSzs=(30 60 90)
maxColSzs=(100 200)
do_rates=(0 0.1)
simThreshs=(0.5)


# backSPIN
   
   # Levels:
   # Level_0_group (always 0)
   # Level_1_group (0 - 1) 
   # Level_2_group (0 - 3)
   # Level_3_group (0 - 7)
   # Level_4_group (0 - 15)

numLevels=(3 4 5 6)
stop_const=(0.5 0.825 1.15)
low_thrs=(0.1 0.2 0.3)



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
                                MEM_caclust=10G
                            fi

                        else
                            MEM_caclust=$MEMORY
                        fi


                        algorithm="caclust_leiden"
                        SCRIPT="${scripts_path}/${algorithm}.R"

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"

                        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --dims $d    \\
   --prune $prune \\
   --NNs $n \\
   --resolution $res \\
   --sim TRUE \\
   --truth $truth \\
   --graph_select $gs \\
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh


                        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh
                                   
                        algorithm="caclust_igraph"
                        SCRIPT="${scripts_path}/${algorithm}.R"

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"


                        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

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
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh


                        algorithm="caclust_spectral"
                        SCRIPT="${scripts_path}/${algorithm}.R"

                        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"
                        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

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
   --usegap TRUE \\
   --graph_select $gs \\
   --nclust NULL \\
   --SNN_mode $snn \\
   --gcKNN $gc \\
   --overlap $ov \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
                               --threads=$THREADS \
                               --memory=$MEM_caclust \
                               -t $MINUTES \
                               bash $tmp_sh


                        counter=$((counter+1))
                        
                        done
                    done
                done                              
            done
        done
    done


#     # end caclust

	#################
	#### divbiclust #
	#################

    algorithm="divbiclust"
    SCRIPT="${scripts_path}/${algorithm}.R"
    
    for maxdiff in ${maxdiffs[@]}; do
        for seedColSz in ${seedColSzs[@]}; do
            for maxColSz in ${maxColSzs[@]}; do
                for simThresh in ${simThreshs[@]}; do
                    for do_rate in ${do_rates[@]};do
    
    
                nm="${algorithm}_${filename}_ntop-${nt}_maxdiff-${maxdiff}_seedColSz-${seedColSz}_maxColSz-${maxColSz}_simThresh-${simThresh}_doRate-${do_rate}"

                tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"
cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $filename \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
    --maxdiff $maxdiff \\
    --seedColSz $seedColSz \\
    --maxColSz $maxColSz \\
    --simThresh $simThresh \\
    --dorate $do_rate \\
    --isdivbiclust TRUE \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                       --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
                       --threads=$THREADS \
                       --memory=$MEMORY \
                       -t $MINUTES \
                       bash $tmp_sh

    	            done
            	done
        	done
    	done
    done


    #########
    # QUBIC #
    #########

    algorithm="QUBIC"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for r in ${r_param[@]}; do
        for q in ${q_param[@]}; do
            for c in ${c_param[@]}; do

                if (($r >= 6))
                then
                    MEM_QUBIC=20G
                else
                    MEM_QUBIC=$MEMORY
                fi

                nm="${algorithm}_${filename}_ntop-${nt}_rparam-${r}_qparam-${q}_cparam-${c}"
                tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_QUBIC
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --r_param $r \\
   --q_param $q \\
   --c_param $c \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

                mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
                       --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
                       --threads=$THREADS \
                       --memory=$MEM_QUBIC \
                       -t $MINUTES \
                       bash $tmp_sh

            done
        done
    done

#     # end QUBIC


    ########
    # s4vd #
    ########

    algorithm="s4vd"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for pv in ${pcerv[@]}; do
        for pu in ${pceru[@]}; do
            for ss in ${ss_thr_min[@]}; do
                for sa in ${ss_thr_add[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_pcerv-${pv}_pceru-${pu}_ssthrmin-${ss}_ssthradd-${sa}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --pcerv $pv \\
   --pceru $pu \\
   --ss_thr_min $ss \\
   --ss_thr_add $sa \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh

                done
            done
        done
    done

    # end s4vd

    #########
    # Plaid #
    #########

    algorithm="Plaid"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for rr in ${rrelease[@]}; do
        for cr in ${crelease[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_rrelease-${rr}_crelease-${cr}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --rrelease $rr \\
   --crelease $cr \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh


        done
    done

    # end Plaid

    ##########
    # Unibic #
    ##########

    algorithm="Unibic"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for tp in ${t_param[@]}; do
        for qd in ${q_discr[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_tparam-${tp}_qdiscr-${qd}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$((MINUTES*2))
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --t_param $tp \\
   --q_discr $qd \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $((MINUTES*2)) \
               bash $tmp_sh


        done
    done

    # end Unibic

    #########
    # Bimax #
    #########

    algorithm="Bimax"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for mr in ${minr[@]}; do
        for mc in ${minc[@]}; do
            # for mac in ${maxc[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_minr-${mr}_minc-${mc}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT  \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --minr $mr \\
   --minc $mc \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh


            # done
        done
    done

    # end Bimax

    #######
    # CCA #
    #######

    algorithm="CCA"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for de in ${delta[@]}; do
        for al in ${alpha[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_delta-${de}_alpha-${al}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --delta $de \\
   --alpha $al \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh


        done
    done

    # end CCA

    ###########
    # Xmotifs #
    ###########

    algorithm="Xmotifs"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for ns in ${ns_param[@]}; do
        for nd in ${nd_param[@]}; do
            for sd in ${sd_param[@]}; do
                for ax in ${alpha_X[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_nsparam-${ns}_ndparam-${nd}_sdparam-${sd}_alphaX-${ax}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --ns_param $ns \\
   --nd_param $nd \\
   --sd_param $sd \\
   --alphaX $ax \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh


                done
            done
        done
    done

#     # end Xmotifs

    ###############
    # BC_Spectral #
    ###############

    algorithm="BC_Spectral"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for nE in ${numEig[@]}; do
        for mrB in ${minr_BCS[@]}; do
            for mcB in ${minc_BCS[@]}; do
                for wV in ${withinVar[@]}; do
                    # for nbB in ${nbest[@]}; do
                    #     for ncs in ${n_clusters[@]}; do


        # nm="${algorithm}_${dataset}_${filename}_ntop-${nt}_numEig-${nE}_minrBCS-${mrB}_mincBCS-${mcB}_withinVar-${wV}_nbest-${nbB}_nclusters-${ncs}"
        nm="${algorithm}_${filename}_ntop-${nt}_numEig-${nE}_minrBCS-${mrB}_mincBCS-${mcB}_withinVar-${wV}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --numEig $nE \\
   --minr_BCS $mrB \\
   --minc_BCS $mcB \\
   --withinVar $wV \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh

                   # --nbest $nbB \
                   # --nclusters $ncs

                    #     done
                    # done
                done
            done
        done
    done

    # end BC_spectral

    ##################
    # QUBIC2-IRISFGM #
    ##################

    algorithm="IRISFGM"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for qCo in ${qcons[@]}; do
        for qOv in ${qoverlap[@]}; do
            for qCm in ${qcmin[@]}; do

        nm="${algorithm}_${filename}_ntop-${nt}_Qconsistency-${qCo}_Qoverlap-${qOv}_Qcmin-${qCm}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=5G
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --Qconsistency $qCo \\
   --Qoverlap $qOv \\
   --Qcmin $qCm \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               --tmpdir=5G \
               -t $MINUTES \
               bash $tmp_sh


            done
        done
    done

    # end IRISFGM

    ##########
    # Seurat #
    ##########

    algorithm="Seurat"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for d in ${dims_seurat[@]}; do
        for n in ${NNs_seurat[@]}; do
            for r in ${res_seurat[@]}; do
                for mp in ${min_perc[@]}; do
                    for lt in ${logfc_thr[@]}; do
                        for rt in ${return_thr[@]}; do



        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_resolution-${r}_minperc-${mp}_logfcthr-${lt}_return_thr-${rt}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --dims $d \\
   --NNs $n \\
   --resolution $r \\
   --logfc_thr $lt \\
   --min_perc $mp \\
   --rthr $rt \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh

                        done
                    done
                done
            done
        done
    done

    # end Seurat

    ############
    # Monocle3 #
    ############

    algorithm="Monocle3"
    SCRIPT="${scripts_path}/${algorithm}.R"

    for d in ${dims_monocle[@]}; do
            for r in ${res_monocle[@]}; do
                for m in ${reduce_method[@]}; do
                    for k in ${NNs_monocle[@]}; do
                        for n in ${ngene_perg[@]}; do



        nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_redm-${m}_resolution-${r}_ngene_pg-${n}_NNs-${k}"
        tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset \\
   --name $nm \\
   --ntop $nt \\
   --sim TRUE \\
   --truth $truth \\
   --nclust $nclust \\
   --dims $d \\
   --resolution $r \\
   --ngene_pg $n \\
   --NNs $k \\
   --redm $m \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

        mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
               --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
               --threads=$THREADS \
               --memory=$MEMORY \
               -t $MINUTES \
               bash $tmp_sh


                    done 
                done
            done
        done
    done

    # end Monocle3

   ############
   # BackSPIN #
   ############


   BS_MINUTES=600

   algorithm="BackSPIN"
   SCRIPT="${scripts_path}/${algorithm}.R"

   for l in ${numLevels[@]}; do
       for s in ${stop_const[@]}; do
           for t in ${low_thrs[@]}; do



       nm="${algorithm}_${filename}_ntop-${nt}_numLevels-${l}_stopconst-${s}_lowthrs-${t}"
       tmp_sh="${here_dir}/bench_sim_${dataset}_${nm}.sh"

cat << EOF > $tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$BS_MINUTES
# END_MXQ

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
  --outdir $OUTDIR  \\
  --file $f \\
  --dataset $dataset \\
  --name $nm \\
  --ntop $nt \\
  --sim TRUE \\
  --truth $truth \\
  --nclust $nclust \\
  --numLevels $l \\
  --stop_const $s \\
  --low_thrs $t \\
&& mv $tmp_sh $here_dir/.done/
EOF
chmod +x $tmp_sh

       mxqsub --stdout="${logdir}/bench_sim_${dataset}_${nm}.stdout.log" \
              --group-name="bench_sim_${dataset}_${filename}_${algorithm}" \
              --threads=$THREADS \
              --memory=$MEMORY \
              -t $BS_MINUTES \
              bash $tmp_sh

           done
       done
   done

   # end BackSPIN



    done
done
fi
