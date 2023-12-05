## Evaluating the robustness of CAbiNet over the choice of number of dimensions

Scripts in this folder runs CAbiNet, Seurat and Monocle3 with number of dimensions of CA/PCA space ranging from 2 to 200 on simulated data sets. The clustering reseluts are then compared to investigate the robustness of algorithms with the dimensionality.

- 1. 'cluster_submit.sh' runs the algirhtms.
  2. Excurating 'find_failed_tasks.sh'. the failed runs will be run again by allocating more running time and memory.
- 3. The output results can then be collated by running 'comm_sh.sh' which calls the R script 'collate_results.R'.
- 4. The evaluation is then visualized by '../../Results/SupplFigures.ipynb'