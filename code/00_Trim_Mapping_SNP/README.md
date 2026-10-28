# Read processing, mapping and variant calling

These five shell scripts were run on the institute SLURM cluster and produced the BAM files,
coverage summaries and VCFs that the downstream modules start from. They are **not** needed to
reproduce modules `01`–`09`, which use the files in `data/`.

| Script | Step |
|---|---|
| `1_trim.sh` | Trimmomatic 0.39 PE trimming (adapter file, `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60`) |
| `2_mapping.sh` | bwa-mem2 against the IId-HLJ 11730 assembly, then `sambamba markdup -r` and `samtools flagstat` |
| `3_coverage.sh` | per-position depth summary with `samtools mpileup`, breadth of coverage and mean depth per chromosome |
| `bcftools-ploidy1.sh` | `bcftools mpileup -a FORMAT/AD,FORMAT/DP` → `bcftools call --ploidy 1 -mv` (haploid samples) |
| `bcftools-ploidy2.sh` | the same with `--ploidy 2` (diploid samples) |

## Before running

The paths below are specific to the cluster on which the pipeline was executed and must be edited:

* `1_trim.sh` — `java -jar .../trimmomatic-0.39.jar`, `ILLUMINACLIP:.../TruSeq3-PE.fa`
* `2_mapping.sh` — bwa-mem2 index prefix `/biodata/.../11730`
* `bcftools-ploidy*.sh` — reference FASTA `/biodata/.../11730_spades_r.fasta`

Each script reads a sample list from the working directory (`SampleID_cut1`, `refID`,
`Sample_path`), one sample (or one BAM path) per line, and writes into sub-directories created by the
scripts (`cleandata/`, `markdup_ref_HLJ/`, `coverage/`, `../vcf_hap/`).

The `#SBATCH` headers request the original partitions, node names and thread counts; delete them when
running the commands outside SLURM.

## Variant set used downstream

The 1,391 markers analysed in the manuscript are the intersection of GATK4 HaplotypeCaller and
BCFtools calls passing `QUAL ≥ 30`, `FORMAT/DP ≥ 10` and `INFO/MQ ≥ 40` (1,138 SNPs + 253 InDels).
The bcftools side of that call is produced by `bcftools-ploidy1.sh` / `bcftools-ploidy2.sh`; the
merged, filtered VCFs are provided in `data/`.
