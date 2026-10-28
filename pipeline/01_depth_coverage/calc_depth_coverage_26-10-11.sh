#!/bin/bash
set -euo pipefail

# ============================================================
# Calculate depth, genome coverage and mapping rate from BAM
# ============================================================

# -------- User configuration --------
REF_FASTA="/share/home/hty/work/IId-cross-Single-Oocyst-26-06-03/IIdA20G1-Harbin-11730-TGS-24-10-01/11730_spades_r.fasta"

BAM_DIR="markdup_ref_HLJ"
OUT_DIR="coverage_stats"

THREADS=8

# 与后续变异分析保持一致
MIN_MAPQ=0
MIN_BASEQ=0

# ------------------------------------

mkdir -p "$OUT_DIR"

# Reference index
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Creating reference index..."
    samtools faidx "$REF_FASTA"
fi

GENOME_SIZE=$(awk '{sum+=$2} END {print sum}' "${REF_FASTA}.fai")

echo "Reference genome: $REF_FASTA"
echo "Reference genome size: $GENOME_SIZE bp"
echo "BAM directory: $BAM_DIR"

OUT_FILE="$OUT_DIR/all_coverage_depth_mapping.txt"

printf "Sample\tGenome_bp\tPrimary_total\tPrimary_mapped\tMapping_pct\tProperly_paired\tProperly_paired_pct\tMean_depth\tSD_depth\tDepth_CV\tCov_1x_bp\tCov_1x_pct\tCov_3x_bp\tCov_3x_pct\tCov_10x_bp\tCov_10x_pct\tCov_30x_bp\tCov_30x_pct\n" \
    > "$OUT_FILE"

shopt -s nullglob
BAMS=( "$BAM_DIR"/*_markdup.bam )

if [ "${#BAMS[@]}" -eq 0 ]; then
    echo "No BAM files found in $BAM_DIR"
    exit 1
fi

for BAM in "${BAMS[@]}"; do

    sample=$(basename "$BAM" _markdup.bam)

    echo "Processing $sample ..."

    # Check BAM
    samtools quickcheck -v "$BAM"

    # Create BAM index if missing
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # --------------------------------------------------------
    # 0. Locate and parse flagstat log
    #    与 BAM 同目录，命名形如 <sample>_stat.log
    #    例如: markdup_ref_HLJ/W585_stat.log
    # --------------------------------------------------------

    STAT_LOG="$(dirname "$BAM")/${sample}_stat.log"

    if [ ! -f "$STAT_LOG" ]; then
        echo "WARNING: flagstat log not found for $sample : $STAT_LOG"
        PP_COUNT="NA"
        PP_PCT="NA"
    else
        # 从形如下面的行中提取：
        # 4691206 + 0 properly paired (10.28% : N/A)
        read -r PP_COUNT PP_PCT < <(
            awk '
                /properly paired/ {
                    count = $1
                    pct = $0
                    sub(/^.*\(/, "", pct)
                    sub(/%.*$/, "", pct)
                    print count, pct
                    exit
                }
            ' "$STAT_LOG"
        )

        # 若未匹配到则设为 NA
        PP_COUNT="${PP_COUNT:-NA}"
        PP_PCT="${PP_PCT:-NA}"
    fi

    # --------------------------------------------------------
    # 1. Mapping statistics
    #
    # 2304 = 256 secondary + 2048 supplementary
    # 2308 = 4 unmapped + 256 secondary + 2048 supplementary
    # --------------------------------------------------------

    PRIMARY_TOTAL=$(samtools view -@ "$THREADS" -c -F 2304 "$BAM")

    PRIMARY_MAPPED=$(samtools view -@ "$THREADS" -c -F 2308 "$BAM")

    MAPPING_PCT=$(awk \
        -v mapped="$PRIMARY_MAPPED" \
        -v total="$PRIMARY_TOTAL" \
        'BEGIN {
            if (total > 0) {
                printf "%.4f", mapped / total * 100
            } else {
                printf "NA"
            }
        }')

    # --------------------------------------------------------
    # 2. Depth and coverage
    #
    # -aa: output every position, including zero-depth bases
    # -Q: minimum mapping quality
    # -q: minimum base quality
    # -d 0: unlimited maximum depth
    # --------------------------------------------------------

    samtools depth \
        -aa \
        -Q "$MIN_MAPQ" \
        -q "$MIN_BASEQ" \
        -d 0 \
        "$BAM" 2>/dev/null |
    awk \
        -v genome="$GENOME_SIZE" \
        -v id="$sample" \
        -v primary_total="$PRIMARY_TOTAL" \
        -v primary_mapped="$PRIMARY_MAPPED" \
        -v mapping_pct="$MAPPING_PCT" \
        -v properly_paired="$PP_COUNT" \
        -v properly_paired_pct="$PP_PCT" '
    BEGIN {
        n=0
        sum=0
        sumsq=0
        cov1=0
        cov3=0
        cov10=0
        cov30=0
    }

    {
        depth=$3

        n++
        sum += depth
        sumsq += depth * depth

        if (depth >= 1)  cov1++
        if (depth >= 3)  cov3++
        if (depth >= 10) cov10++
        if (depth >= 30) cov30++
    }

    END {
        if (n == 0) {
            print id "\tNA"
            exit
        }

        mean_depth = sum / n
        sd_depth = sqrt((sumsq / n) - (mean_depth * mean_depth))

        if (mean_depth > 0) {
            depth_cv = sd_depth / mean_depth
        } else {
            depth_cv = 0
        }

        cov1_pct  = cov1  / genome * 100
        cov3_pct  = cov3  / genome * 100
        cov10_pct = cov10 / genome * 100
        cov30_pct = cov30 / genome * 100

        printf "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",
            id,
            genome,
            primary_total,
            primary_mapped,
            mapping_pct,
            properly_paired,
            properly_paired_pct,
            mean_depth,
            sd_depth,
            depth_cv,
            cov1,
            cov1_pct,
            cov3,
            cov3_pct,
            cov10,
            cov10_pct,
            cov30,
            cov30_pct
    }' >> "$OUT_FILE"

done

echo
echo "Done."
echo "Results saved to:"
echo "$OUT_FILE"
echo
column -t -s $'\t' "$OUT_FILE" | head -n 10