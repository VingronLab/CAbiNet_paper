
Die Daten sind simuliert mit Splatter (https://bioconductor.org/packages/release/bioc/html/splatter.html) und  library size, mean gene expression, etc wurden geschätzt auf der Basis der 'Zeisel Brain Data' (https://www.science.org/doi/10.1126/science.aaa1934).

Ich hab 1000 Zellen, 10000 Gene und 6 Cluster (Proportionen: 0.25, 0.1, 0.1, 0.2, 0.3, 0.05) simuliert.

Preprocessing etc. ist für die logcounts auch bereits gemacht, aber raw counts sind auch enthalten.

Die 6 Dateien hier im Ordner unterscheiden sich nur an der stärke der DE in den clustern.

dePROB -> Wahrscheinlichkeit dass ein Gen im Cluster differential exprimiert ist.

Wenn ein Gen differential exprimiert ist wird es multipliziert mit einem Faktor der aus einer log-normalen Verteilung kommt:
defacLOC -> Mittelwert der log-normal Verteilung
defSCALE -> Varianz der log-normal Verteilung

Die Dateien kannst du laden via readRDS() und sind SingleCellExperiment Container.

Die Metainformationen für die Gene ist enthalten in rowData(simulated_data), für die Zellen in colData(simulated_data).

