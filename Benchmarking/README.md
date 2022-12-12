
# Benchmarking

This folder contains:

- Folder "algorithms" with R scripts for running the (bi-)clustering algorithms mentioned in the manuscript.
- The R script "setup.R", which prepares the data and calls the script for the needed algorithm.
- The bash scripts "cluster_submit_real.sh" and "cluster_submit_sim.sh" for experimental and simulated data respectively.
- The R script "collate_results.R", which collects the output of the benchmarking and saves it in a single table.

The bash scripts for cluster submission contain the parameters used for the benchmarking.
As we cannot anticipate the cluster environment you will run the scripts on, the exact commands will likely need to be changed. 
Currently the "cluster_submit_(real/sim).sh" scripts use the inhouse version of qsub, mxqsub.

"setup.R" gets called from the cluster submission scripts with the appropriate parameters, loads the data and calls the script for the (bi-)clustering algorithm.
As the algorithms are run in parallel on a cluster, "collate_results.R" is needed to combine the output .csv files to a single table.

