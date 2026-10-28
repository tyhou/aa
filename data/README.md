# Data

Input tables, variant calls and intermediate files for the analyses in `code/`. The directory layout
mirrors the paths the scripts expect under `/data`, so mounting this folder at `/data` is enough to
run the pipeline.

Everything here is released under CC0 1.0 (see `LICENSE` in this folder).

## Contents

| Path | What it is | Used by |
|---|---|---|
| `01_Het_Fws/fws/*_hap_overlap_1391.vcf` | 1,391-marker haploid VCFs for the expansion lines and for the F1 single oocysts/sporozoites | `Fws-hap-*.R` (converted to GDS on the fly) |
| `01_Het_Fws/het/*` | PLINK `--het` and `--hardy` output and the diploid VCF it was computed from | `het-cov.R` |
| `01_Het_Fws/het-cov.txt` | Heterozygosity and breadth of coverage at ≥1×, ≥10× and ≥30× for every sample (source of the heterogeneity/coverage panel) | `het-cov.R` |
| `02_PCA/{F1SOSS,SOLine}/` | PLINK `--pca` eigenvectors/eigenvalues and sample group files | PCA scripts |
| `02_PCA/PCA_QC_metadata.xlsx` | Depth, mapping rate, heterogeneity, group, DPI and batch for each sample — the covariates tested against the PCs | `SOLine_overlap_1391_PC_QC_heatmap.R` |
| `03_Variant_distribution/intersection_{snps,indels}_positions.txt` | Positions of the 1,138 SNPs and 253 InDels shared by GATK4 and BCFtools | `variant_distribution_plot.R` |
| `03_Variant_distribution/Marker_1065.txt` | Positions of the earlier 1,065-SNP set, for comparison | — |
| `04_Allele_Frequency/F1SOSS/`, `VAF_Density/` | Haploid/diploid VCFs for the F1 bulk progeny, single oocysts and single sporozoites | VAF and density scripts |
| `04_Allele_Frequency/SOLine/` | Allele-frequency table of the expansion lines, the haploid VCF and the per-line heterogeneity input | crossover and haplotype-stability scripts |
| `05_Genetic_Map/` | Crossover events (absolute VAF difference > 0.8), LDhat `Rho.stat.final.txt`, chromosome sizes and the crossover–heterozygosity table | genetic-map scripts |
| `06_Fst_Pi/` | Windowed F<sub>ST</sub> and Pi (1,000-bp windows, 100-bp steps) and the IId 11730 GFF3 | `ROD-fst-plot-win1000-step100.R` |
| `07_Line3_8_11_24_Allele_Frequency/` | Allele frequencies and the list of compared samples | Mann–Whitney script |
| `08_XP-EHH/phasing/` | Phased per-chromosome VCFs for all samples and for the high- and low-virulence groups | `xpehh_manhattan.R` |
| `09_Upset_Plot/list*.txt` | The four candidate-gene sets | `Upset.R` |
| `11_AlphaFold/` | AlphaFold 3 prediction archive for CpDHR1-HLJ: full-data JSON (per-residue pLDDT, PAE) and model CIF | `11_AlphaFold/*.R` |

## Marker sets

Two marker sets appear in this folder. The 1,391-marker set (1,138 SNPs + 253 InDels, GATK4 ∩
BCFtools) is the one used for the analyses reported in the revised manuscript; the 1,065-SNP set is
the earlier call set, kept because the PLINK heterozygosity and PCA commands in `01_Het_Fws/het/`
and `02_PCA/pca_dip.sh` were written against it. The PLINK outputs computed from the 1,391-marker
VCFs are provided directly in `02_PCA/*/` (`*_dip_overlap_1391_PCA.eigenvec`), so the PCA scripts
run without repeating that step.

## Not included

* raw FASTQ/BAM files (NCBI BioProject `TODO-BIOPROJECT-ACCESSION`);
* the reference assembly (IIdA20G1-HLJ 11730);
* GROMACS trajectories for the molecular-dynamics analysis (756 MB; available from the Code Ocean
  capsule).
