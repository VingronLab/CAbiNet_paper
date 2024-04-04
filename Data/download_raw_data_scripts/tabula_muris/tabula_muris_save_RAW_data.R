
library(TabulaMurisData)

sce <- TabulaMurisSmartSeq2()
sce <- sce[, sce$tissue == "Limb_Muscle"]
sce <- sce[, !is.na(sce$cell_ontology_class)]

saveRDS(sce, "../../discussed_data/raw/tabula_muris_RAW.rds")
