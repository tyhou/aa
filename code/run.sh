#!/usr/bin/env bash
# Master script: runs every analysis module in dependency order.
#
# Container (recommended):
#   docker run --rm -v "$PWD/data":/data -v "$PWD/code":/code -v "$PWD/results":/results \
#     cpar-revision bash /code/run.sh
#
# Native R session:
#   BASE_CODE="$PWD/code" RESULTS_BASE="$PWD/results" bash code/run.sh
#   (after rewriting the /data and /results prefixes, see README "Option B")

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_CODE="${BASE_CODE:-$SCRIPT_DIR}"
RESULTS_BASE="${RESULTS_BASE:-/results}"

# shared helper functions for the allele-frequency scripts
export ALLELEFREQ_CODE_DIR="${ALLELEFREQ_CODE_DIR:-$BASE_CODE/04_Allele_Frequency}"

if [ ! -d /data ]; then
    echo "WARNING: /data does not exist. Mount data/ at /data, or rewrite the input paths (README, Option B)." >&2
fi

scripts=(
    "01_Het_Fws/het_coverage_plot/het-cov.R"
    "01_Het_Fws/Het_dotplot.R"
    "01_Het_Fws/fws/Fws-hap-SOLine-1391.R"
    "01_Het_Fws/fws/Fws-hap-SOSS-1391.R"
    "02_PCA/F1SOSS_PCA_overlap_1391.R"
    "02_PCA/SOLine_PCA_overlap_1391.R"
    "02_PCA/SOLine_overlap_1391_PC_QC_heatmap.R"
    "03_Variant_distribution/variant_distribution_plot.R"
    "04_Allele_Frequency/VAF-F1-bulk-SOSS-density_overlap_1391.R"
    "04_Allele_Frequency/F1SOSS_alleleFreq_1391.R"
    "04_Allele_Frequency/SOLine_alleleFreq_crossover_1391_Het_lastDPI.R"
    "04_Allele_Frequency/SOLine_timecourse_haplotype_stability.R"
    "05_Genetic_Map/crossover_het_r2.R"
    "05_Genetic_Map/CrossoverCounts.R"
    "05_Genetic_Map/recombinationMap.R"
    "06_Fst_Pi/ROD-fst-plot-win1000-step100.R"
    "07_Line3_8_11_24_Allele_Frequency/alleleFreq_Line3_8_11_vs_24.R"
    "08_XP-EHH/xpehh_manhattan.R"
    "09_Upset_Plot/Upset.R"
    "10_Drift_Concordance/SOLine_parallel_dp.R"
    "11_AlphaFold/cgd5_4090_hlj-pLDDT.R"
    "11_AlphaFold/cgd5_4090_hlj-pae.R"
    "11_AlphaFold/cgd5_4090_hlj-pLDDT-legend.R"
    "12_Sensitivity_W689del/ROD-fst-plot-fws0.5-W689del.R"
)

total=${#scripts[@]}
index=0

for script_rel in "${scripts[@]}"; do
    index=$((index + 1))
    script="$BASE_CODE/$script_rel"
    result_dir="$RESULTS_BASE/$(dirname "$script_rel")"

    mkdir -p "$result_dir"

    echo "=== [$index/$total] $script_rel"
    start=$SECONDS
    Rscript "$script" "$@"
    echo "    done in $((SECONDS - start)) s -> $result_dir"
done

echo "=== all $total analyses finished; outputs are under $RESULTS_BASE"
