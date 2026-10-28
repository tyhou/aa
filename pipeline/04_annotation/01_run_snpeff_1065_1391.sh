#!/usr/bin/env bash
set -euo pipefail

BCFTOOLS="/share/home/hty/software/bcftools-1.10/bcftools"

SNPEFF_DIR="/share/home/hty/software/snpEff"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPEFF_CONFIG="${SNPEFF_DIR}/snpEff.config"

SNPEFF_DATABASE="11730_TGS"

REFERENCE="11730_spades_r.fasta"

OLD_VCF="Parent_hap_filter_maf0.05_geno0.02.vcf"
NEW_VCF="intersection.vcf.gz"

ROOT_DIR="$(pwd)"
OUTPUT_DIR="${ROOT_DIR}/snpeff"

NORMALIZED_DIR="${OUTPUT_DIR}/normalized"
OLD_OUTPUT_DIR="${OUTPUT_DIR}/old_1065"
NEW_OUTPUT_DIR="${OUTPUT_DIR}/new_1391"

mkdir -p \
    "${NORMALIZED_DIR}" \
    "${OLD_OUTPUT_DIR}" \
    "${NEW_OUTPUT_DIR}"

echo "========================================"
echo "Checking programs and input files"
echo "========================================"

if [[ ! -x "${BCFTOOLS}" ]]; then
    echo "ERROR: bcftools cannot be executed:"
    echo "${BCFTOOLS}"
    exit 1
fi

if [[ ! -s "${SNPEFF_JAR}" ]]; then
    echo "ERROR: SnpEff jar was not found:"
    echo "${SNPEFF_JAR}"
    exit 1
fi

if [[ ! -s "${SNPEFF_CONFIG}" ]]; then
    echo "ERROR: SnpEff configuration was not found:"
    echo "${SNPEFF_CONFIG}"
    exit 1
fi

for file in \
    "${REFERENCE}" \
    "${REFERENCE}.fai" \
    "${OLD_VCF}" \
    "${NEW_VCF}"
do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: missing or empty input file: ${file}"
        exit 1
    fi
done

command -v java >/dev/null 2>&1 || {
    echo "ERROR: Java was not found."
    exit 1
}

echo "Java:"
java -version 2>&1 | head -n 3

echo
echo "SnpEff:"
java -jar "${SNPEFF_JAR}" \
    -version 2>&1 |
    head -n 5

echo
echo "SnpEff database: ${SNPEFF_DATABASE}"

echo
echo "========================================"
echo "Normalizing the two VCF files"
echo "========================================"

# Split multiallelic variants and verify REF against the reference genome.
"${BCFTOOLS}" norm \
    -f "${REFERENCE}" \
    -m -any \
    -Oz \
    -o "${NORMALIZED_DIR}/old_1065.norm.vcf.gz" \
    "${OLD_VCF}"

"${BCFTOOLS}" norm \
    -f "${REFERENCE}" \
    -m -any \
    -Oz \
    -o "${NORMALIZED_DIR}/new_1391.norm.vcf.gz" \
    "${NEW_VCF}"

"${BCFTOOLS}" index \
    -f \
    -t \
    "${NORMALIZED_DIR}/old_1065.norm.vcf.gz"

"${BCFTOOLS}" index \
    -f \
    -t \
    "${NORMALIZED_DIR}/new_1391.norm.vcf.gz"

OLD_NORMALIZED_COUNT=$(
    "${BCFTOOLS}" view \
        -H \
        "${NORMALIZED_DIR}/old_1065.norm.vcf.gz" |
    wc -l |
    tr -d ' '
)

NEW_NORMALIZED_COUNT=$(
    "${BCFTOOLS}" view \
        -H \
        "${NORMALIZED_DIR}/new_1391.norm.vcf.gz" |
    wc -l |
    tr -d ' '
)

echo "Normalized old variants: ${OLD_NORMALIZED_COUNT}"
echo "Normalized new variants: ${NEW_NORMALIZED_COUNT}"

echo
echo "========================================"
echo "Annotating old 1065 SNP set"
echo "========================================"

(
    cd "${OLD_OUTPUT_DIR}"

    java \
        -Xmx8g \
        -jar "${SNPEFF_JAR}" \
        ann \
        -c "${SNPEFF_CONFIG}" \
        -v \
        -no-upstream \
        -no-downstream \
        -stats "snpEff_summary.html" \
        "${SNPEFF_DATABASE}" \
        "${NORMALIZED_DIR}/old_1065.norm.vcf.gz" \
        > "old_1065.snpeff.vcf" \
        2> "old_1065.snpeff.log"
)

echo
echo "========================================"
echo "Annotating new 1391 SNP+InDel set"
echo "========================================"

(
    cd "${NEW_OUTPUT_DIR}"

    java \
        -Xmx8g \
        -jar "${SNPEFF_JAR}" \
        ann \
        -c "${SNPEFF_CONFIG}" \
        -v \
        -no-upstream \
        -no-downstream \
        -stats "snpEff_summary.html" \
        "${SNPEFF_DATABASE}" \
        "${NORMALIZED_DIR}/new_1391.norm.vcf.gz" \
        > "new_1391.snpeff.vcf" \
        2> "new_1391.snpeff.log"
)

echo
echo "========================================"
echo "Compressing and indexing annotated VCFs"
echo "========================================"

"${BCFTOOLS}" view \
    -Oz \
    -o "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf.gz" \
    "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf"

"${BCFTOOLS}" view \
    -Oz \
    -o "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf.gz" \
    "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf"

"${BCFTOOLS}" index \
    -f \
    -t \
    "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf.gz"

"${BCFTOOLS}" index \
    -f \
    -t \
    "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf.gz"

echo
echo "========================================"
echo "Checking annotated records"
echo "========================================"

OLD_ANNOTATED_COUNT=$(
    "${BCFTOOLS}" view \
        -H \
        "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf.gz" |
    wc -l |
    tr -d ' '
)

NEW_ANNOTATED_COUNT=$(
    "${BCFTOOLS}" view \
        -H \
        "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf.gz" |
    wc -l |
    tr -d ' '
)

OLD_WITH_ANN=$(
    "${BCFTOOLS}" query \
        -f '%INFO/ANN\n' \
        "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf.gz" |
    awk '$0 != "." {count++} END {print count+0}'
)

NEW_WITH_ANN=$(
    "${BCFTOOLS}" query \
        -f '%INFO/ANN\n' \
        "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf.gz" |
    awk '$0 != "." {count++} END {print count+0}'
)

echo "Old input records:      ${OLD_NORMALIZED_COUNT}"
echo "Old annotated records:  ${OLD_ANNOTATED_COUNT}"
echo "Old records with ANN:   ${OLD_WITH_ANN}"

echo
echo "New input records:      ${NEW_NORMALIZED_COUNT}"
echo "New annotated records:  ${NEW_ANNOTATED_COUNT}"
echo "New records with ANN:   ${NEW_WITH_ANN}"

if [[ "${OLD_NORMALIZED_COUNT}" -ne "${OLD_ANNOTATED_COUNT}" ]]; then
    echo "ERROR: old VCF record count changed during annotation."
    exit 1
fi

if [[ "${NEW_NORMALIZED_COUNT}" -ne "${NEW_ANNOTATED_COUNT}" ]]; then
    echo "ERROR: new VCF record count changed during annotation."
    exit 1
fi

echo
echo "========================================"
echo "SnpEff annotation completed"
echo "========================================"

echo "Old annotated VCF:"
echo "${OLD_OUTPUT_DIR}/old_1065.snpeff.vcf.gz"

echo
echo "New annotated VCF:"
echo "${NEW_OUTPUT_DIR}/new_1391.snpeff.vcf.gz"

echo
echo "Old SnpEff report:"
echo "${OLD_OUTPUT_DIR}/snpEff_summary.html"

echo
echo "New SnpEff report:"
echo "${NEW_OUTPUT_DIR}/snpEff_summary.html"