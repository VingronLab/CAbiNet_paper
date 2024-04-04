
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc

input_mat = pd.read_table("../../discussed_data/raw/GSE189981_d50_organoids_cnt_matrix.txt.gz")
metadata = pd.read_table("../../discussed_data/raw/GSE189981_d50_organoids_metadata.txt.gz")


cnt_mat_info = input_mat.loc[:, "cell_id":"protocol"]
cnt_mat = input_mat.loc[:, "MIR1302-2HG":]

in_metadata = input_mat.cell_id.isin(metadata.cell_id)
cnt_mat = cnt_mat[in_metadata]

adata_merged = sc.AnnData(cnt_mat)
adata_merged.obs = metadata

adata_merged.write("../../discussed_data/raw/brain_organoids_adata_merged_RAW.h5ad")

