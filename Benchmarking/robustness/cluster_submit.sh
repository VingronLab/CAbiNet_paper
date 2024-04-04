#!/bin/bash

# If you run this script from the directory its located in
# set working to the parent directory. Otherwise change accordingly.
cd ../

# NOTE: Monocle3 requires specific libraries. Change path to where your libraries are installed.
LD_LIBRARY_VAR=../renv/local/lib:../renv/local/include
GDAL_DATA_VAR=../renv/local/share/gdal

scripts_path="./algorithms/"
outdir=./results
mkdir -p $outdir

OUTDIR="${outdir}/out/"
mkdir -p $OUTDIR

# NOTE: Make sure you ran sim_generation.R first!
# (also from the data_scalability folder!)
# If you prefer to only run either the simulated data with 6 clusters or 20 clusters, comment the respective parts of the code below.

# 6 clusters
indir_a='../Data/sim_data/raw/data_scalability/'
# datasets_=`ls -d ${indir}/*`

# 20 clusters
indir_b='../Data/sim_data/raw/robustness/'
# datasets=$(ls -d ${indir}/pbmc3k*deprob-0.1_defacLoc-0.75_defacScale-0.75*)

# Run 6 and 20 cluster datasets at the same time. Comment if running only one of them!
datasets=$(ls -d ${indir_a}/*/ ${indir_b}/pbmc3k*deprob-0.1_defacLoc-0.75_defacScale-0.75*/)

for dataset in ${datasets[@]}; do

	dataset_bn=$(basename $dataset)

	logdir="${outdir}/log/${dataset_bn}"
	here_dir="${outdir}/sh/${dataset_bn}"
	mkdir -p $logdir
	mkdir -p $here_dir
	mkdir -p $here_dir/.done/

	files="${dataset}/*.rds"

	ntop=(500) #500 for datasets with 2000 genes
	nclust=6
	truth='Group'

	####################################################################################

	# CAclust - 108
	res=1
	dims=(2 4 6 8 10 20 30 40 50 100 150 200)
	NNs=(30)
	overlap=(0.1)
	graph_select=(FALSE)
	SNN_mode=("all")
	calc_gckNN=(FALSE)
	prune=$(echo 'print(1/15)' | python3)

	# Seurat
	res_seurat=(1)
	dims_seurat=(2 4 6 8 10 20 30 40 50 100 150 200)
	NNs_seurat=(30)
	min_perc=(0.25)
	logfc_thr=(0.25)
	return_thr=(0.05)

	# Monocle 108
	res_monocle=(1)
	dims_monocle=(2 4 6 8 10 20 30 40 50 100 150 200)
	NNs_monocle=(30)
	reduce_method=("PCA")
	ngene_perg=(50)

	for f in ${files[@]}; do

		filename=$(basename $f .rds)

		# Extract the substring after "ncell-"
		substring=$(echo $filename | grep -o "ncell[^_]*")
		# Extract the number using regular expression
		number=$(echo "$substring" | grep -o '[0-9]*')

		if [[ "${substring}" == "ncell-1e+05" ]]; then
			number=100000
		fi

		# Set runtime, memory etc. depending on the number.
		if [[ $number -le 10000 ]]; then
			THREADS=6
			MEMORY=100G
			MINUTES=600
			CA_MEM=200G
		elif [[ $number -gt 10000 && $number -le 60000 ]]; then
			THREADS=12
			MEMORY=300G
			MINUTES=1080
			CA_MEM=600G
		else
			THREADS=32
			MEMORY=600G
			MINUTES=1440
			CA_MEM=600G
		fi

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
											MEM_caclust=$CA_MEM
										fi

									else
										MEM_caclust=$MEMORY
									fi

									echo $MEM_caclust

									algorithm="caclust_igraph"
									SCRIPT="${scripts_path}/${algorithm}.R"

									nm="${algorithm}_${filename}_ntop-${nt}_dims-${d}_NNs-${n}_gs-${gs}_SNN-${snn}_gcKNN-${gc}_overlap-${ov}"

									tmp_sh="${here_dir}/bench_sim_${nm}.sh"
									cat <<EOF >$tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEM_caclust
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset_bn \\
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

									mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
										--group-name="bench_sim_${filename}_${algorithm}" \
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

			#     end caclust
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

									tmp_sh="${here_dir}/bench_sim_${nm}.sh"

									cat <<EOF >$tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset_bn \\
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

									mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
										--group-name="bench_sim_${filename}_${algorithm}" \
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

								tmp_sh="${here_dir}/bench_sim_${nm}.sh"

								cat <<EOF >$tmp_sh
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# t=$MINUTES
# END_MXQ

export LD_LIBRARY_PATH=$LD_LIBRARY_VAR:\$LD_LIBRARY_PATH
export GDAL_DATA=$GDAL_DATA_VAR:\$GDAL_DATA

trap 'echo ERROR_TIMEOUT >&2' SIGXCPU

Rscript-4.2.1 $SCRIPT   \\
   --outdir $OUTDIR  \\
   --file $f \\
   --dataset $dataset_bn \\
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

								mxqsub --stdout="${logdir}/bench_sim_${nm}.stdout.log" \
									--group-name="bench_sim_${filename}_${algorithm}" \
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

			#     # end Monocle3

		done
	done
done
