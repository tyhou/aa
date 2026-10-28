# Single-oocyst genomics of a *Cryptosporidium parvum* IId-HLJ × IId-GD genetic cross

Code, input data and computational environment for the bioinformatic analyses of the study
*"High recombination and intra-oocyst heterogeneity shape virulence in Cryptosporidium parvum, with
CpDHR1 defined as a novel virulence determinant"*.

The virulent IIdA20G1-HLJ strain (tdTomato) and the avirulent IIdA19G1-GD strain (mNeonGreen) were
crossed in IFN-γ knockout mice. Dual-fluorescent F1 progeny oocysts were sorted by flow cytometry and
sequenced individually, both as whole oocysts and as single sporozoites after whole-genome
amplification, and 15 single-oocyst expansion lines were established. The scripts below reproduce the
genomic analyses of that material: intra-oocyst heterogeneity, the recombination map, the
allele-frequency contrast between high- and low-virulence lines, and the intersection of candidate
gene sets that nominated CpDHR1 (cgd5_4090).

---

## Repository layout

```
code/          analysis scripts, one directory per analysis module
pipeline/      cluster shell/Python steps that produced the tables in data/
data/          input tables, variant calls and intermediate files used by the scripts
environment/   Dockerfile, postInstall and package versions exported from Code Ocean
```

`data/` mirrors the paths the scripts expect under `/data`, and the scripts write to `/results`;
both are mounted into the container when the pipeline is run (see [Quick start](#quick-start)).

---

## Requirements

**R 4.5.1** with the packages pinned in `environment/Dockerfile`: `SeqArray`, `SeqVarTools`,
`GenomicRanges`, `moimix`, `rehh`, `vcfR`, `tidyverse`, `patchwork`, `ggforce`, `ggpubr`, `UpSetR`,
`pheatmap`, `scatterplot3d`, `openxlsx`, `writexl`, `viridis`, `R.utils` and their dependencies. The
exact versions are fixed in the Dockerfile (`remotes::install_version()`) and in `postInstall`.

---

## Quick start

```bash
docker build -t cpar-revision environment/
mkdir -p results
docker run --rm \
  -v "$PWD/data":/data \
  -v "$PWD/code":/code \
  -v "$PWD/results":/results \
  cpar-revision bash /code/run.sh
```




