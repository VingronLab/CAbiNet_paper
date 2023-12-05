
# Preprocessing

This folder contains:

- The R scripts to process any data set (`data_preprocessing.R` and `data_preprocessing_splatter.R` for real and simulated data respectively).
- The R script to pre-process the brain organoid data set (`data_preprocessing_brain_organoids.R`), which is necessary as the brain organoids need less stringent filtering of mitochondrial reads.
- The bash scripts to submit the preprocessing scripts to a computing cluster (`submit_to_cluster_preprocessing.sh` and `submit_to cluster_preprocessing_simdata.sh`).
- For the brain organoid and tabula muris data set we provide the the `run_preprocessing_brain_organoids.sh` and `run_preprocessing_tabula_muris.sh` files that ensure pre-processing exactly as reported in the manuscript.
- The directory `sim_data_generation`, which will generate the simulated data.

The bash scripts for cluster submission contain the parameters used for the analysis presented in the manuscript.
As we cannot anticipate the cluster environment you will run the scripts on, the exact commands will likely need to be changed. 
Currently the `submit_to_cluster_preprocessing(\_simdata).sh` scripts use the inhouse version of qsub, mxqsub.

# Generation of simulated data

As described in the manuscript, the simulated data is generated based on either the `zeisel` or `pbmc3k` datasets. These files are automatically downloaded when running the `make_splatter_data.R` script.
After generation of the simulated data, it should be pre-processed similarly to the experimental data with the `submit_to_cluster_preprocessing_simdata.sh` script.
