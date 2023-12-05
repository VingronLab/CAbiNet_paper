
# CAbiNet - Data sets used in the manuscript

All the raw/preprocessed datasets, the benchmarking results and data sets discussed in paper can be found in the online data repository https://zenodo.org/record/7433294#.Y6HA5rTMI11.

- The scripts to download the brain organoid and tabula muris data sets can be found in the `download_raw_data_scripts` folder.
- Experimental data sets, both raw and preprocessed, should be deposited into the respective folder in `real_data`.
- Simulated data sets should be copied either into the `sim_data/raw` or `sim_data/preprocessed` folders. The datases in folder `sim_data/raw/data_scalability` are for the evaluation of the scalability and robustness of CAbiNet. The codes for the evaluation can be found from `Benchmarking/data_scalability` and `Benchmarking/robustness`.

In the folder `real_data/`, we provided 9 data sets we downloaded from the public resources as listed in the Table 2 in our manuscript and in the supplementary `data_source.xlsx` file. Some of the data sets are subsetted and assigned expert annotated cell types to the metadata following the codes in the file:
- `DataFormatting.ipynb`

Tips for data downloading:
- Simply download the data from the online repository and copy the unzipped files into the ./Data directory. This will provide you with all the data to reproduce our results. For reproducing the results through the Jupyter notebooks provided in `Results`, you only need the `*_data/preprocessed` files.
- If you would like to also run your own pre-processing, you can remove the `real_data/preprocessed`, `sim_data/preprocessed` as well as the `discussed_data/preprocessed` files, as they will be regenerated during the pipeline. Make sure there are files in the `*_data/raw_data` subdirectories.
