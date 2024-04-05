
# CAbiNet - code accompanying the manuscript

This code reproduces the figures in our manuscript "CAbiNet: Joint visualization of cells and genes based on a gene-cell graph" that can be found here: **Link to paper**.


We split the code into several sections to make it easier to reproduce the findings:

* Preprocessing: 
    - Simulated data sets generated with Splatter from two scRNA-seq data sets.
    - Preprocessing of all data sets: Normalization, dimension reduction, feature selection, etc. ...
* Benchmarking:
    - Benchmarking of simulated data sets.
    - Benchmarking of experimental scRNA-seq data sets.
    - Scripts for collating the results.
    - Evaluation of the scalability and robustness of CAbiNet.
* Results:
    - Jupyter-notebooks containing the results shown in figures (Fig.2-6) in our manuscript as well as the supplementary figures.
* Data:
    - Scripts that download and standardize the raw data discussed in Results.
* SupplementaryMaterial:
    - Interactive plots.


If you wish to redo the whole analysis, you can simply download all the data sets from the repository https://zenodo.org/records/10260709 and https://zenodo.org/records/10932001 to the empty folder `Data`, and follows the 3 sections (Preprocessing -> Benchmarking -> Results) mentioned above. If you simply extract the downloaded data into the `./Data` directory you will have all the data necessary to reproduce our results. The version of "CAbiNet" that used to reproduce results in our paper can be found from https://figshare.com/articles/software/CAbiNet_Joint_clustering_and_visualization_of_cells_and_genes_for_single-cell_transcriptomics/23276402, the file named as "CAbiNet-main.zip".

In case you would like to only use parts of the data, here are some tips:

- Download the RAW data into the correct folder under either `Data/sim_data/raw`, `Data/real_data/raw` or `Data/discussed_data/raw`.
- Alternatively you can reproduce the results starting with the preprocessed data and follow the scripts in `Benchmarking` folder. Simply download preprocessed data into `Data/sim_data/preprocessed`, `Data/real_data/preprocessed` or `Data/discussed_data/preprocessed`

The easiest way to reproduce the findings in the manuscript is to download the archive from the online repository and to unpack all the data into the `Data` directory and to then run the jupyter notebooks in the `Results` folder.
Necessary data files to reproduce plots from the manuscript:
* Data/discussed_data/preprocessed 
* Data/sim_data/preprocessed
* Data/real_data/preprocessed
* Benchmarking/results (precomputed benchmarking results)

Some scripts may need to be adapted to your computational environment, as they were created with the server structure of the Max Planck Institute for Molecular Genetics in mind.
