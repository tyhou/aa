# Computational environment

`Dockerfile` and `postInstall` are exported verbatim from the Code Ocean capsule used for the
revised manuscript. Together they define the exact system on which every number and figure in the
manuscript was produced:

* base: `codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04`
  (Ubuntu 22.04, R 4.5.1, RStudio Server 2025.05.1);
* R packages installed at fixed versions with `remotes::install_version()`
  (`tidyverse 2.0.0`, `vcfR 1.16.0`, `rehh 3.2.3`, `UpSetR 1.4.1`, `patchwork 1.3.2`,
  `pheatmap 1.0.13`, `scatterplot3d 0.3-45`, `openxlsx 4.2.9`, `writexl 2.0.1`, `viridis 0.6.5`,
  `R.utils 2.13.0`, and the `tidyverse` dependency set);
* Bioconductor 3.22 packages (`SeqArray`, `SeqVarTools`, `GenomicRanges`, `IRanges`, `S4Vectors`,
  `BiocParallel`, `BiocStyle`) installed by `postInstall`;
* `moimix` (Fws) from `https://plasmogenepi.r-universe.dev`;
* `ggchicklet` from a fixed GitHub commit.
