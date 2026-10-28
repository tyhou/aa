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

**Command-line tools** (module `00`): Trimmomatic 0.39, bwa-mem2, samtools, sambamba, bcftools,
PLINK 1.9 (`01`–`02`). The container installs the R side; the aligners are expected on `PATH` if the
raw reads are re-processed from FASTQ.

---

## Quick start

### Option A — container (recommended)

```bash
docker build -t cpar-revision environment/
mkdir -p results
docker run --rm \
  -v "$PWD/data":/data \
  -v "$PWD/code":/code \
  -v "$PWD/results":/results \
  cpar-revision bash /code/run.sh
```

`run.sh` runs the 24 R scripts in dependency order and writes each module's output to
`results/<module>/`.

### Option B — local R session

The scripts use the absolute paths `/data` and `/results` so that they run unchanged inside the
container. To run them natively, rewrite those two prefixes once and then execute the master script:

```bash
grep -rl '"/data\|"/results' code --include='*.R' | \
  xargs sed -i "s|\"/data|\"$PWD/data|g; s|\"/results|\"$PWD/results|g"
mkdir -p results
BASE_CODE="$PWD/code" RESULTS_BASE="$PWD/results" bash code/run.sh
```

All analyses are pure R and run on a laptop; only module `00` needs a cluster-scale machine.

---

## Data

**Included** (~33 MB): the 1,391-marker haploid and diploid VCFs, the SNP and InDel position lists,
PLINK `.het`/`.eigenvec` outputs, windowed F<sub>ST</sub> and Pi tables, 48 phased per-chromosome
VCFs, the IId 11730 GFF3 annotation, the LDhat `Rho.stat.final.txt`, the crossover tables, and the
AlphaFold 3 prediction archive for CpDHR1-HLJ (full-data JSON and model CIF, 15 MB). These are
sufficient to re-run modules `01`–`12` without re-processing any reads.

**Not included**:

* raw FASTQ and BAM files — deposited in NCBI (BioProject `TODO-BIOPROJECT-ACCESSION`);
* the reference assembly (IIdA20G1-HLJ 11730, TGS assembly) — public accession
  `TODO-REFERENCE-ACCESSION`;
* the molecular-dynamics trajectories (756 MB of GROMACS output; not part of this repository).

### The 1,391-marker variant set

The marker set is the intersection of GATK4 HaplotypeCaller and BCFtools calls passing
`QUAL ≥ 30`, `FORMAT/DP ≥ 10` and `INFO/MQ ≥ 40`, and comprises 1,138 SNPs and 253 InDels; 1,017 of
the 1,065 SNPs from the earlier set are retained. The filtered VCFs and the two position lists used
by module `03` are in `data/`. Upstream bcftools calls can be regenerated with
`code/00_Trim_Mapping_SNP/bcftools-ploidy1.sh` (haploid) and `bcftools-ploidy2.sh` (diploid).


