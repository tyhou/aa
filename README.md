# Single-oocyst genomics of a *Cryptosporidium parvum* IId-HLJ × IId-GD genetic cross

Code, input data and computational environment for the bioinformatic analyses of the study
*"High recombination and intra-oocyst heterogeneity shape virulence in Cryptosporidium parvum, with
CpDHR1 defined as a secondary virulence determinant"*.

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

## Analysis modules

| Module | Analysis | Inputs (in `data/`) | Main outputs |
|---|---|---|---|
| `00_Trim_Mapping_SNP` | Adapter/quality trimming (Trimmomatic 0.39), mapping to the IId-HLJ 11730 assembly (bwa-mem2), duplicate marking (sambamba), per-base depth summaries (samtools), variant calling (bcftools, `--ploidy 1` and `--ploidy 2`) | raw FASTQ (not included), reference FASTA | BAM files, `all_coverage_depth_final.txt`, `all_hap.vcf.gz` |
| `01_Het_Fws` | Marker heterozygosity (PLINK `--het`) and Fws per sample (SeqArray + moimix), and the heterozygosity–coverage relationship | `01_Het_Fws/**` (1,391-marker VCFs, PLINK `.het`) | `Fws_*.txt`, `heterozygosity_coverage_plot.pdf/png` |
| `02_PCA` | PCA of the expansion lines and of the F1 single oocysts/sporozoites (PLINK `--pca`), group- and heterogeneity-coloured projections, scree plots and a PC–covariate heatmap | `02_PCA/**` | `*_PCA.eigenvec`, PCA and QC figures |
| `03_Variant_distribution` | Density of the 1,138 SNPs and 253 InDels across the eight chromosomes | `03_Variant_distribution/intersection_*_positions.txt` | variant-density figures |
| `04_Allele_Frequency` | Alternative-allele frequency profiles of the parents, the F1 bulk progeny, single oocysts and single sporozoites; crossover calls between adjacent markers (absolute VAF difference > 0.8); early- versus late-sample breakpoint stability; VAF density distributions | `04_Allele_Frequency/**` | allele-frequency tables, `crossover.xlsx`, VAF and density figures |
| `05_Genetic_Map` | Physical-versus-genetic map, crossover counts, telomere-distance distribution and the crossover–heterozygosity correlation | `05_Genetic_Map/**` | `map-15-Lines.csv`, crossover statistics, figures |
| `06_Fst_Pi` | Windowed F<sub>ST</sub> and nucleotide diversity (1,000-bp windows, 100-bp steps) between the high- and low-virulence groups, with outlier windows | `06_Fst_Pi/**` | `merge_pi_Fst-1000win-100step.xlsx`, F<sub>ST</sub>/Pi figures |
| `07_Line3_8_11_24_Allele_Frequency` | Mann–Whitney contrast of allele frequencies between lines 3, 8, 11 and line 24 | `07_Line3_8_11_24_Allele_Frequency/**` | significant-SNP tables |
| `08_XP-EHH` | XP-EHH between the high- and low-virulence groups from the phased per-chromosome VCFs (rehh), Manhattan plots | `08_XP-EHH/**` | `*.xpehh.txt`, Manhattan figures |
| `09_Upset_Plot` | Intersection of the four candidate-gene sets (highly polymorphic, XP-EHH, F<sub>ST</sub>/Pi, allele-frequency diversity) | `09_Upset_Plot/list*.txt` | `gene_lists_intersection.xlsx`, UpSet plot |
| `10_Drift_Concordance` | Direction concordance of allele-frequency change between the first and last sample of each line (same fraction), binned, with gene annotation of the concordant variants | `04_Allele_Frequency/SOLine/*.tsv`, `06_Fst_Pi/IId_11730_genome.gff3` | `parallel_dp_summary.csv`, `same_frac_bin_pct.pdf`, `candidate_snp_genes.csv` |
| `11_AlphaFold` | pLDDT profile from the AlphaFold 3 CIF/JSON output, PAE matrix and heatmap, domain-wise confidence, pLDDT colour legend | `11_AlphaFold/fold_*_hlj_full_data_0.json`, `*_model_0.cif` | pLDDT and PAE figures, `pae_matrix.csv`, domain statistics |
| `12_Sensitivity_W689del` | Sensitivity re-plot of the F<sub>ST</sub>/Pi outlier analysis with the W689 line excluded | `06_Fst_Pi/*.windowed.*`, GFF3 | ROD/F<sub>ST</sub> density figures, merged table |

---

## Upstream command-line steps

`pipeline/` holds the cluster shell and Python steps that produced the tables in `data/`:

* depth, breadth of coverage at ≥1×, ≥3×, ≥10× and ≥30×, and mapping statistics per BAM
  (`01_depth_coverage`);
* Kraken2 composition of the reads that do not map to the parasite genome (`02_kraken2`);
* per-sample splitting of the merged call set and extraction of the 1,065-SNP and 1,391-marker
  subsets, plus PLINK heterozygosity (`03_variant_selection`);
* snpEff annotation, the SNP-set Venn diagram and the highly polymorphic gene set (`04_annotation`);
* the LDhat recombination pipeline and the conversion of rho to cM (`05_ldhat`);
* windowed F<sub>ST</sub>, Pi and Tajima's D (`06_fst_pi`);
* sequenceLDhot, hotspot extraction and motif parsing (`07_sequenceldhot`).

They are not needed to run `code/` — the outputs they produce are already in `data/` — but they
document exactly how those inputs were generated. See `pipeline/README.md` for the script-to-output
map and the software versions.

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


