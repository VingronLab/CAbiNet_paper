
# Preprocessing

This folder contains:
- the R scripts to process any data set ("data_preprocessing.R" and "data_preprocessing_splatter.R" for real and simulated data respectively)
- the bash scripts to submit the preprocessing scripts to a computing cluster ("submit_to_cluster_preprocessing.sh" and "submit_to cluster_preprocessing_simdata.sh").
- empty directory ("ExperimentalData" and "sim_data") in which the input data should be placed.
- The directory "sim_data_generation", which will generate the simulated data.

The bash scripts for cluster submission contain the parameters used for the analysis presented in the manuscript.
As we cannot anticipate the cluster environment you will run the scripts on, the exact commands will likely need to be changed. 
Currently the "submit_to_cluster_preprocessing(\_simdata).sh" scripts use the inhouse version of qsub, mxqsub.

# Generation of simulated data

As described in the manuscript, the simulated data is generated based on either the "zeisel" or "pbmc3k" datasets. Instructions on where to download the files can be found in the manuscript or in the README.md in the root folder of the repository.
The zeisel and pbmc3k SingleCellExperiment objects should be placed as .rds files into the "./sim_data/data/preprocessed/" directory.

After generation of the simulated data, it should be pre-processed similarly to the experimental data with the "submit_to_cluster_preprocessing_simdata.sh" script.
