# Upstream command-line steps

These are the command-line steps that turn the sequencing reads into the tables in `data/`. They ran
on the institute SLURM cluster; the paths (`/share/home/hty`, `/biodata/hty`) and the `#SBATCH`
headers are cluster-specific and must be edited before re-use.

The R analyses in `code/` start from the outputs of this folder, so the two are complementary: this
folder documents how the inputs were produced, `code/` reproduces the figures and statistics.

| Folder | Script | What it does | Output used by |
|---|---|---|---|
| `01_depth_coverage` | `calc_depth_coverage_26-10-11.sh` | Per-BAM mapping statistics and breadth of coverage at ≥1×, ≥3×, ≥10× and ≥30× (`samtools depth`, flagstat), plus mean depth and depth CV | `data/01_Het_Fws/het-cov.txt`, Fig. 2a, Supplementary Tables 1–3 |
| `02_kraken2` | `parse_kraken_reports_100pct.py` | Collapses Kraken2 reports into a mutually exclusive 100% composition (Unclassified, Eukaryota, Bacteria, Archaea, Viruses, Artificial, other) with *Cryptosporidium* as a non-additive sub-annotation | contamination assessment, Supplementary Table 2 |
| `03_variant_selection` | `01_filter_SO_hap_260708.sh`, `01_filter_SO_dip_260708.sh` | Splits the merged call set per sample (`bcftools view -S`) and extracts the 1,065-SNP and 1,391-marker subsets (`bcftools view -R`) | all downstream modules |
| `03_variant_selection` | `02_het_calculation.sh` | PLINK `--het`/`--hardy` for both marker sets (F1 bulk, single oocysts/sporozoites, expansion lines) | `data/01_Het_Fws/het/`, heterogeneity values |
| `04_annotation` | `01_run_snpeff_1065_1391.sh` | snpEff annotation of the old (1,065) and new (1,391) variant sets against the `11730_TGS` database | Response Table 2 |
| `04_annotation` | `02_extract_snpeff_annotations.py` | Pulls the coding consequences out of the annotated VCFs | Response Table 2 |
| `04_annotation` | `plot_SNP_venn_1065_vs_1391.R` | Venn diagram of the two SNP sets | Response Fig. 3 |
| `04_annotation` | `filter_high_polymorphism_genes.R` | Genes ranked by coding-variant density (highly polymorphic set) | `data/09_Upset_Plot/list4`, Table 1, Fig. 4g |
| `05_ldhat` | `run_ldhat_pipeline_nseq30.sh` | vcftools LDhat format → `lkgen -nseq 30` → `interval` (1,000,000 iterations) → `stat`, merged to `Rho.stat.final.txt` | `data/05_Genetic_Map/Rho.stat.final.txt` |
| `05_ldhat` | `convert_Rho2cM.sh` | Converts LDhat rho to genetic distance in cM with *N*<sub>e</sub> = 39.604 | kb/cM estimate |
| `05_ldhat` | `add_span_and_summary_1391.py`, `add_span_and_summary_1065.py` | Adds interval spans and per-chromosome summaries to the rho tables | genetic-map tables |
| `05_ldhat` | `plot_compare_cumulative.py` | Compares cumulative genetic maps (needs pandas/matplotlib) | map comparison figure |
| `06_fst_pi` | `fst-pi-tajimaD-by-pop-1000win-100step.sh` | `vcftools` windowed F<sub>ST</sub> (`--weir-fst-pop`), Pi and Tajima's D, 1,000-bp windows / 100-bp steps | `data/06_Fst_Pi/*`, Fig. 4f |
| `07_sequenceldhot` | `sequenceLDhot_hty_26-10-10.sh` | Runs sequenceLDhot on the phased per-chromosome data and summarises hotspots | `data/05_Genetic_Map` hotspots, Supplementary Fig. 10 |
| `07_sequenceldhot` | `run_hotspot_subseq.sh` | `seqkit subseq` extracts each hotspot window from the assembly | motif input |
| `07_sequenceldhot` | `extract_genes_in_hotspot.py` | Maps hotspots to genes | hotspot gene list |
| `07_sequenceldhot` | `meme-motif-extract.py` | Parses MEME output for the enriched motifs | Supplementary Fig. 10 |

## Software

The scripts call, besides the R stack in `environment/`: bcftools 1.10, bgzip/tabix, vcftools,
PLINK 1.9, samtools, snpEff (database `11730_TGS`), Kraken2 (reports are the input here), LDhat,
sequenceLDhot, seqkit, MEME. Python scripts use the standard library only, except
`plot_compare_cumulative.py` (pandas, numpy, matplotlib).

## What is not scripted here

The step that merged the GATK4 and BCFtools calls into `all_hap.vcf.gz` / `all_dip.vcf.gz` and then
applied `QUAL ≥ 30`, `FORMAT/DP ≥ 10`, `INFO/MQ ≥ 40` to produce the `*.id.vcf.gz` files was run
interactively and is not preserved as a script. The filtering rule is stated in the Methods, and the
resulting marker sets are provided as `data/03_Variant_distribution/intersection_snps_positions.txt`
and `intersection_indels_positions.txt`, so the two calls can be reproduced and checked position by
position.
