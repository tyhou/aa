#!/usr/bin/env python3

import csv
import gzip
import os
from collections import Counter, defaultdict


INPUT_VCFS = {
    "Old_1065": (
        "snpeff/old_1065/"
        "old_1065.snpeff.vcf.gz"
    ),
    "New_1391": (
        "snpeff/new_1391/"
        "new_1391.snpeff.vcf.gz"
    )
}

OUTPUT_DIR = "snpeff/tables"

ANN_COLUMNS = [
    "Allele",
    "Effect",
    "Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS_c",
    "HGVS_p",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Distance",
    "Warnings"
]

IMPACT_SCORE = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0
}

# Used only when several annotations have the same IMPACT.
EFFECT_SCORE = {
    "chromosome_number_variation": 100,
    "exon_loss_variant": 95,
    "frameshift_variant": 90,
    "stop_gained": 88,
    "stop_lost": 86,
    "start_lost": 84,
    "splice_acceptor_variant": 82,
    "splice_donor_variant": 82,
    "rare_amino_acid_variant": 80,
    "missense_variant": 70,
    "inframe_insertion": 68,
    "inframe_deletion": 68,
    "protein_altering_variant": 65,
    "splice_region_variant": 60,
    "synonymous_variant": 50,
    "stop_retained_variant": 48,
    "start_retained_variant": 48,
    "5_prime_UTR_variant": 40,
    "3_prime_UTR_variant": 40,
    "intron_variant": 30,
    "non_coding_transcript_variant": 25,
    "upstream_gene_variant": 20,
    "downstream_gene_variant": 20,
    "intergenic_region": 10
}


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8"
        )

    return open(
        path,
        "r",
        encoding="utf-8"
    )


def parse_info(info_text):
    result = {}

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = True

    return result


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    if len(ref) < len(alt):
        return "Insertion"

    if len(ref) > len(alt):
        return "Deletion"

    if len(ref) == len(alt) and len(ref) > 1:
        return "MNV"

    return "Complex"


def functional_region(effect):
    effects = effect.split("&")

    if "intergenic_region" in effects:
        return "Intergenic"

    if any(
        "UTR" in item
        for item in effects
    ):
        return "UTR"

    if "intron_variant" in effects:
        return "Intron"

    if any(
        item in effects
        for item in [
            "splice_acceptor_variant",
            "splice_donor_variant",
            "splice_region_variant"
        ]
    ):
        return "Splice_region"

    coding_effects = {
        "synonymous_variant",
        "missense_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "start_retained_variant",
        "stop_retained_variant",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
        "coding_sequence_variant"
    }

    if any(
        item in coding_effects
        for item in effects
    ):
        return "Coding"

    if any(
        item in effects
        for item in [
            "upstream_gene_variant",
            "downstream_gene_variant"
        ]
    ):
        return "Upstream_downstream"

    return "Other_gene_related"


def annotation_score(row):
    impact = IMPACT_SCORE.get(
        row["Impact"],
        0
    )

    effect_scores = [
        EFFECT_SCORE.get(effect, 0)
        for effect in row["Effect"].split("&")
    ]

    effect_score = max(
        effect_scores,
        default=0
    )

    return impact, effect_score


def parse_vcf(dataset, vcf_file):
    rows = []

    with open_text(vcf_file) as handle:
        for line in handle:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0]
            position = fields[1]
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filter_value = fields[6]
            info_text = fields[7]

            info = parse_info(info_text)
            ann_text = info.get("ANN", "")

            if not ann_text:
                rows.append({
                    "Dataset": dataset,
                    "CHROM": chrom,
                    "POS": position,
                    "ID": variant_id,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": quality,
                    "FILTER": filter_value,
                    "Variant_Type": variant_type(
                        ref, alt
                    ),
                    "Annotation_Number": 0,
                    **{
                        column: ""
                        for column in ANN_COLUMNS
                    },
                    "Functional_Region": (
                        "Unannotated"
                    )
                })

                continue

            for annotation_number, ann in enumerate(
                ann_text.split(","),
                start=1
            ):
                values = ann.split("|")

                if len(values) < len(ANN_COLUMNS):
                    values.extend(
                        [""] * (
                            len(ANN_COLUMNS) -
                            len(values)
                        )
                    )

                if len(values) > len(ANN_COLUMNS):
                    values = values[
                        :len(ANN_COLUMNS)
                    ]

                ann_values = dict(
                    zip(
                        ANN_COLUMNS,
                        values
                    )
                )

                row = {
                    "Dataset": dataset,
                    "CHROM": chrom,
                    "POS": position,
                    "ID": variant_id,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": quality,
                    "FILTER": filter_value,
                    "Variant_Type": variant_type(
                        ref, alt
                    ),
                    "Annotation_Number": (
                        annotation_number
                    )
                }

                row.update(ann_values)

                row["Functional_Region"] = (
                    functional_region(
                        row["Effect"]
                    )
                )

                rows.append(row)

    return rows


def write_table(path, rows, fieldnames):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def write_counter_summary(
    path,
    counter,
    category_name
):
    totals = defaultdict(int)

    for (dataset, category), count in counter.items():
        totals[dataset] += count

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t"
        )

        writer.writerow([
            "Dataset",
            category_name,
            "Count",
            "Percentage"
        ])

        for dataset, category in sorted(counter):
            count = counter[
                (dataset, category)
            ]

            total = totals[dataset]

            percentage = (
                100 * count / total
                if total > 0
                else 0
            )

            writer.writerow([
                dataset,
                category,
                count,
                f"{percentage:.4f}"
            ])


def main():
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True
    )

    all_rows = []

    for dataset, vcf_file in INPUT_VCFS.items():
        if not os.path.isfile(vcf_file):
            raise FileNotFoundError(
                f"Input VCF not found: "
                f"{vcf_file}"
            )

        dataset_rows = parse_vcf(
            dataset,
            vcf_file
        )

        all_rows.extend(dataset_rows)

        print(
            f"{dataset}: "
            f"{len(dataset_rows)} ANN records"
        )

    base_fields = [
        "Dataset",
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "Variant_Type",
        "Annotation_Number"
    ]

    all_fields = (
        base_fields +
        ANN_COLUMNS +
        ["Functional_Region"]
    )

    all_effects_file = os.path.join(
        OUTPUT_DIR,
        "snpeff_all_effects.tsv"
    )

    write_table(
        all_effects_file,
        all_rows,
        all_fields
    )

    # Select one representative annotation per variant.
    grouped_rows = defaultdict(list)

    for row in all_rows:
        key = (
            row["Dataset"],
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"]
        )

        grouped_rows[key].append(row)

    primary_rows = []

    for key, rows in grouped_rows.items():
        primary_row = max(
            rows,
            key=annotation_score
        )

        primary_rows.append(primary_row)

    primary_rows.sort(
        key=lambda row: (
            row["Dataset"],
            row["CHROM"],
            int(row["POS"])
        )
    )

    primary_effects_file = os.path.join(
        OUTPUT_DIR,
        "snpeff_primary_effects.tsv"
    )

    write_table(
        primary_effects_file,
        primary_rows,
        all_fields
    )

    effect_counts = Counter(
        (
            row["Dataset"],
            row["Effect"] or "Unannotated"
        )
        for row in primary_rows
    )

    impact_counts = Counter(
        (
            row["Dataset"],
            row["Impact"] or "Unannotated"
        )
        for row in primary_rows
    )

    region_counts = Counter(
        (
            row["Dataset"],
            row["Functional_Region"]
        )
        for row in primary_rows
    )

    variant_type_counts = Counter(
        (
            row["Dataset"],
            row["Variant_Type"]
        )
        for row in primary_rows
    )

    write_counter_summary(
        os.path.join(
            OUTPUT_DIR,
            "snpeff_effect_summary.tsv"
        ),
        effect_counts,
        "Effect"
    )

    write_counter_summary(
        os.path.join(
            OUTPUT_DIR,
            "snpeff_impact_summary.tsv"
        ),
        impact_counts,
        "Impact"
    )

    write_counter_summary(
        os.path.join(
            OUTPUT_DIR,
            "snpeff_region_summary.tsv"
        ),
        region_counts,
        "Functional_Region"
    )

    write_counter_summary(
        os.path.join(
            OUTPUT_DIR,
            "snpeff_variant_type_summary.tsv"
        ),
        variant_type_counts,
        "Variant_Type"
    )

    gene_counts = Counter(
        (
            row["Dataset"],
            row["Gene_ID"],
            row["Gene_Name"]
        )
        for row in primary_rows
        if row["Gene_ID"] not in ("", ".")
    )

    gene_summary_file = os.path.join(
        OUTPUT_DIR,
        "snpeff_gene_summary.tsv"
    )

    with open(
        gene_summary_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t"
        )

        writer.writerow([
            "Dataset",
            "Gene_ID",
            "Gene_Name",
            "Variant_Count"
        ])

        for key in sorted(gene_counts):
            writer.writerow([
                key[0],
                key[1],
                key[2],
                gene_counts[key]
            ])

    print()
    print("Written:")
    print(all_effects_file)
    print(primary_effects_file)
    print(gene_summary_file)


if __name__ == "__main__":
    main()