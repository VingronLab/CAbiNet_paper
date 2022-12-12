
# CAbiNet - code accompanying the manuscript

This code reproduces the figures in our manuscript "CAbiNet:Joint visualization of cells and genes based on a gene-cell association graph" that can be found here: ADD LINK TO PAPER.


We split the code into 3 sections to make it easier to reproduce the findings:

* Preprocessing: 
    - Simulated data sets generated with Splatter from two scRNA-seq data sets.
    - Preprocessing of all data sets: Normalization, dimension reduction, feature selection, etc. ...
* Benchmarking:
    - Benchmarking of simulated data sets.
    - Benchmarking of experimental scRNA-seq data sets.
    - Scripts for collecting results.
* Results:
    - Jupyter-notebooks containing the results shown in figures (Fig.2-6) in our manuscript. 

If you wish to redo the whole analysis, you can download all the raw data sets from the repository (LINK HERE!!!!) to the empty folder Data/rawdata, and follows 3 sections (Preprocessing -> Benchmarking -> Results) mentioned above. Tips: download rawdata into
* Data/: 
    - rawdata/  (remeber to unzip the zip-files in folder rawdata/simulated/)

Or you can reproduce the results starting with the preprocessed data and follow the scripts in Benchmarking/ folder. Tips: download preprocessed data into
* Data/: 
    - preprocessed/ (remeber to unzip the zip-files in folder preprocessed/simulated/)

The easiest way to reproduce the findings in the manuscript is to download the benchmarking/ and disscussed_data folders into directory Data/ and then run the jupyter notebooks in the "Results" folder.
* Data/: 
    - benchmarking/  
    - discussed_data/

Some scripts may need to be adapted to your computational environment, as they were created with the server structure of the Max Planck Institute for Molecular Genetics in mind.