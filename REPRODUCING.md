# Reproducing the analyses

This repository is the GitHub mirror of the Code Ocean capsule used for the revised manuscript. The
container reproduces the capsule environment (R 4.5.1 and pinned package versions), and `code/run.sh`
runs the 24 analysis scripts in dependency order — the 18 capsule scripts plus the drift-concordance,
AlphaFold and sensitivity modules added for the GitHub release.

## 1. Build the environment

```bash
docker build -t cpar-revision environment/
```

The build installs R 4.5.1 and every package at the version recorded in
`environment/Dockerfile`. It takes 15–30 min on a normal connection and needs ~4 GB of disk.

## 2. Run the pipeline

```bash
mkdir -p results
docker run --rm \
  -v "$PWD/data":/data \
  -v "$PWD/code":/code \
  -v "$PWD/results":/results \
  cpar-revision bash /code/run.sh
```

The scripts read from `/data` and write to `/results`, which is why those two directories are
mounted at the container root. Each module writes to `results/<module>/`, mirroring the layout of
`code/`.

A full run takes about 10–20 min on a laptop; the Fws and XP-EHH steps dominate.

## 3. Expected outputs

| Directory | Files |
|---|---|
| `results/01_Het_Fws/` | `fws/Fws_*.txt`, `het_coverage_plot/heterozygosity_coverage_plot.pdf` |
| `results/02_PCA/` | PCA projections, scree plots, `PC_QC_heatmap` |
| `results/03_Variant_distribution/` | SNP and InDel density figures |
| `results/04_Allele_Frequency/` | allele-frequency tables, crossover tables, VAF density figures |
| `results/05_Genetic_Map/` | genetic map, crossover counts, telomere-distance figures |
| `results/06_Fst_Pi/` | `merge_pi_Fst-1000win-100step.xlsx`, F<sub>ST</sub>/Pi figures |
| `results/07_Line3_8_11_24_Allele_Frequency/` | significant-SNP tables |
| `results/08_XP-EHH/` | XP-EHH tables and Manhattan plots |
| `results/09_Upset_Plot/` | `gene_lists_intersection.xlsx`, UpSet plot |
| `results/10_Drift_Concordance/` | `parallel_dp_summary.csv`, `same_frac_bin_pct.pdf`, `candidate_snp_genes.csv` |
| `results/11_AlphaFold/` | pLDDT profile and by-domain table, `pae_matrix.csv`, PAE heatmap, colour legend |
| `results/12_Sensitivity_W689del/` | F<sub>ST</sub>/Pi re-plot excluding line W689, merged table |

## 4. Running without a container

See "Option B" in the README: rewrite the `/data` and `/results` prefixes in the scripts with the
provided `sed` command, then run `bash code/run.sh`. R 4.5.1 and the packages listed in
`environment/Dockerfile` must be installed.

## Troubleshooting

* **`registry.codeocean.com` cannot be reached while building.** Use the fallback definition:
  `docker build -t cpar-revision -f environment/Dockerfile.portable environment/`. It builds the same
  R 4.5.1 environment from `rocker/r-ver` instead of the Code Ocean base image.
* **A module cannot find its input.** Confirm that `data/` is mounted at `/data` and that the file
  names are unchanged; the scripts address inputs by relative path under `/data`.
* **`seqVCF2GDS` complains about the `.gds` file.** The GDS files are generated on the fly from the
  VCFs; they are not tracked in the repository and may be deleted at any time.
