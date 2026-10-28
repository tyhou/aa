#!/bin/bash
#SBATCH --job-name=fst-pi-tajimaD
#SBATCH -p sonmi
#SBATCH -n 24
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-0
##SBATCH -mem=200G
##SBATCH -t 72:00:00
ulimit -s unlimited
ulimit -l unlimited

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# ============================================================
# Configurable inputs
#
# Default:
#   sbatch fst-pi-tajimaD-by-pop-1000win-100step.sh
#
# Use another VCF without editing this script:
#   sbatch fst-pi-tajimaD-by-pop-1000win-100step.sh another_file.vcf
#   sbatch fst-pi-tajimaD-by-pop-1000win-100step.sh another_file.vcf.gz
# ============================================================
VCF_FILE="${1:-SOLine_dip_overlap_1391.vcf}"
HIGH_POP_FILE="${HIGH_POP_FILE:-highpop.txt}"
LOW_POP_FILE="${LOW_POP_FILE:-lowpop.txt}"
WINDOW_SIZE="${WINDOW_SIZE:-1000}"
WINDOW_STEP="${WINDOW_STEP:-100}"

# Derive an output prefix automatically from the VCF filename.
VCF_NAME=$(basename "$VCF_FILE")
VCF_STEM=${VCF_NAME%.vcf.gz}
VCF_STEM=${VCF_STEM%.vcf}
OUTPUT_DIR="${OUTPUT_DIR:-${VCF_STEM}_fst_pi_tajimaD}"
mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$VCF_FILE" ]]; then
    echo "ERROR: VCF file not found or empty: $VCF_FILE" >&2
    exit 1
fi
if [[ ! -s "$HIGH_POP_FILE" ]]; then
    echo "ERROR: population file not found or empty: $HIGH_POP_FILE" >&2
    exit 1
fi
if [[ ! -s "$LOW_POP_FILE" ]]; then
    echo "ERROR: population file not found or empty: $LOW_POP_FILE" >&2
    exit 1
fi
if ! command -v vcftools >/dev/null 2>&1; then
    echo "ERROR: vcftools is not available in PATH" >&2
    exit 1
fi

# Automatically choose the correct vcftools input option.
if [[ "$VCF_FILE" == *.vcf.gz ]]; then
    VCF_INPUT=(--gzvcf "$VCF_FILE")
else
    VCF_INPUT=(--vcf "$VCF_FILE")
fi

echo "Input VCF : $VCF_FILE"
echo "High group: $HIGH_POP_FILE"
echo "Low group : $LOW_POP_FILE"
echo "Window    : ${WINDOW_SIZE} bp; step ${WINDOW_STEP} bp"
echo "Output dir: $OUTPUT_DIR"

# Windowed Weir-Fst.
vcftools "${VCF_INPUT[@]}" \
    --weir-fst-pop "$HIGH_POP_FILE" \
    --weir-fst-pop "$LOW_POP_FILE" \
    --fst-window-size "$WINDOW_SIZE" \
    --fst-window-step "$WINDOW_STEP" \
    --out "$OUTPUT_DIR/${VCF_STEM}_fst_${WINDOW_SIZE}win_${WINDOW_STEP}step"

# Per-site Weir-Fst.
vcftools "${VCF_INPUT[@]}" \
    --weir-fst-pop "$HIGH_POP_FILE" \
    --weir-fst-pop "$LOW_POP_FILE" \
    --out "$OUTPUT_DIR/${VCF_STEM}_fst_per_site"

# Extract each population, then calculate nucleotide diversity and Tajima's D.
for pop in high low; do
    if [[ "$pop" == "high" ]]; then
        POP_FILE="$HIGH_POP_FILE"
    else
        POP_FILE="$LOW_POP_FILE"
    fi

    KEEP_PREFIX="$OUTPUT_DIR/${VCF_STEM}_keep_${pop}"

    vcftools "${VCF_INPUT[@]}" \
        --keep "$POP_FILE" \
        --recode \
        --recode-INFO-all \
        --out "$KEEP_PREFIX"

    vcftools --vcf "${KEEP_PREFIX}.recode.vcf" \
        --window-pi "$WINDOW_SIZE" \
        --window-pi-step "$WINDOW_STEP" \
        --out "$OUTPUT_DIR/${VCF_STEM}_${pop}_pi_${WINDOW_SIZE}win_${WINDOW_STEP}step"

    vcftools --vcf "${KEEP_PREFIX}.recode.vcf" \
        --TajimaD "$WINDOW_SIZE" \
        --out "$OUTPUT_DIR/${VCF_STEM}_${pop}_TajimaD_${WINDOW_SIZE}win"
done

echo "vcftools analyses finished successfully."
exit
