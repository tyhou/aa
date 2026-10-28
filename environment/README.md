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

## Build

```bash
docker build -t cpar-revision .
```

## Fallback when the Code Ocean registry is unreachable

`Dockerfile.portable` builds the same R 4.5.1 environment from `rocker/r-ver:4.5.1` with the same
package versions, and does not depend on the Code Ocean registry:

```bash
docker build -t cpar-revision -f Dockerfile.portable .
```

Beyond the capsule's package list it adds `jsonlite`, `data.table` and `RColorBrewer`, which the
AlphaFold module (`code/11_AlphaFold`) needs; the capsule image already provided them as dependencies.

Differences from the capsule image: the RStudio Server layer is absent (not needed to run the
scripts), and the package installation is expressed as one `install_version()` call per package
rather than being split between the Dockerfile and `postInstall`. The R version and every package
version are identical. This file was added for the GitHub release and is not part of the original
capsule export; validate it once with a full `code/run.sh` run before citing it.

## Command-line tools

Modules `01`–`09` need R only. Module `00` additionally needs Trimmomatic 0.39, bwa-mem2, samtools,
sambamba and bcftools, which are **not** installed by this image (the reads were processed on the
institute cluster); install them separately if the FASTQ files are to be re-processed.
