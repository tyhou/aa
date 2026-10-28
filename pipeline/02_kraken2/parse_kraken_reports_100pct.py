#!/usr/bin/env python3
"""Create a coarse, mutually exclusive 100% composition from Kraken2 reports.

Cryptosporidium is included in Eukaryota. Its abundance is also reported as a
non-additive annotation, so it must not be added again to the main percentages.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


MAIN_CATEGORIES = [
    "Unclassified",
    "Eukaryota",
    "Bacteria",
    "Archaea",
    "Viruses",
    "Artificial_sequences",
    "Other_classified",
]

TAXIDS = {
    "Unclassified": "0",
    "Eukaryota": "2759",
    "Bacteria": "2",
    "Archaea": "2157",
    "Viruses": "10239",
    "Artificial_sequences": "81077",
    "Cryptosporidium": "5806",
}


def percentages_summing_to_100(counts, total):
    """Round to two decimals with the largest-remainder method (sum=100.00)."""
    if total == 0:
        return {category: 0.0 for category in MAIN_CATEGORIES}
    exact_units = {
        category: counts[category] * 10000 / total for category in MAIN_CATEGORIES
    }
    units = {category: int(exact_units[category]) for category in MAIN_CATEGORIES}
    remaining = 10000 - sum(units.values())
    order = sorted(
        MAIN_CATEGORIES,
        key=lambda category: exact_units[category] - units[category],
        reverse=True,
    )
    for category in order[:remaining]:
        units[category] += 1
    return {category: units[category] / 100 for category in MAIN_CATEGORIES}


def parse_report(report_file):
    """Read clade counts for disjoint top-level taxa and annotations."""
    clade_by_taxid = {}
    species = defaultdict(int)
    root_classified = None

    with report_file.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 6:
                print(f"Warning: skipped malformed line {line_number} in {report_file.name}")
                continue
            _, clade_text, direct_text, rank, taxid, raw_name = parts[:6]
            try:
                clade_reads = int(clade_text)
                direct_reads = int(direct_text)
            except ValueError:
                print(f"Warning: skipped non-numeric line {line_number} in {report_file.name}")
                continue

            name = raw_name.strip()
            clade_by_taxid[taxid] = clade_reads
            if taxid == "1" and rank == "R":
                root_classified = clade_reads
            if rank.startswith("S") and direct_reads > 0:
                species[(taxid, name)] += direct_reads

    if root_classified is None:
        raise ValueError(f"Root row (taxid 1, rank R) not found in {report_file.name}")

    unclassified = clade_by_taxid.get(TAXIDS["Unclassified"], 0)
    total = root_classified + unclassified
    counts = {
        category: clade_by_taxid.get(TAXIDS[category], 0)
        for category in MAIN_CATEGORIES
        if category != "Other_classified"
    }

    selected_classified = sum(
        counts[category]
        for category in MAIN_CATEGORIES
        if category not in {"Unclassified", "Other_classified"}
    )
    counts["Other_classified"] = root_classified - selected_classified
    if counts["Other_classified"] < 0:
        raise ValueError(
            f"Selected top-level clades overlap or exceed root in {report_file.name}."
        )
    if sum(counts.values()) != total:
        raise ValueError(f"Internal count check failed for {report_file.name}")

    crypto_reads = clade_by_taxid.get(TAXIDS["Cryptosporidium"], 0)
    euk_reads = counts["Eukaryota"]
    if crypto_reads > euk_reads:
        raise ValueError(f"Cryptosporidium exceeds Eukaryota in {report_file.name}")
    annotations = {
        "Cryptosporidium_within_Eukaryota_reads": crypto_reads,
        "Cryptosporidium_%_of_total": crypto_reads / total * 100 if total else 0,
        "Cryptosporidium_%_of_Eukaryota": crypto_reads / euk_reads * 100 if euk_reads else 0,
        "Non_Cryptosporidium_Eukaryota_reads": euk_reads - crypto_reads,
        "Non_Cryptosporidium_Eukaryota_%_of_total": (
            (euk_reads - crypto_reads) / total * 100 if total else 0
        ),
    }
    return counts, total, annotations, species


def main():
    parser = argparse.ArgumentParser(
        description="Summarize Kraken2 reports into coarse categories totaling 100%."
    )
    parser.add_argument("--report-dir", default=".", type=Path)
    parser.add_argument("--output-dir", default=".", type=Path)
    args = parser.parse_args()

    reports = sorted(args.report_dir.glob("*.kraken.report"))
    if not reports:
        raise SystemExit(f"No *.kraken.report files found in: {args.report_dir.resolve()}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    top10_rows = []
    for report in reports:
        sample = report.name.removesuffix(".kraken.report")
        print(f"Processing {sample} ...")
        counts, total, annotations, species = parse_report(report)
        percentages = percentages_summing_to_100(counts, total)

        row = {"Sample": sample, "Total_input_sequences": total}
        for category in MAIN_CATEGORIES:
            row[f"{category}_reads"] = counts[category]
            row[f"{category}_%"] = f"{percentages[category]:.2f}"
        row["Main_percentage_sum"] = f"{sum(percentages.values()):.2f}"
        for key, value in annotations.items():
            row[key] = f"{value:.2f}" if isinstance(value, float) else value
        summary_rows.append(row)

        ranked_species = sorted(species.items(), key=lambda item: item[1], reverse=True)
        for rank_number, ((taxid, name), reads) in enumerate(ranked_species[:10], 1):
            top10_rows.append({
                "Sample": sample, "Rank": f"Top{rank_number}", "Species": name,
                "Reads": reads, "TaxID": taxid,
            })

    summary_path = args.output_dir / "contamination_summary_100pct.csv"
    top10_path = args.output_dir / "top10_species_all_samples_100pct.csv"
    with summary_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    with top10_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["Sample", "Rank", "Species", "Reads", "TaxID"]
        )
        writer.writeheader()
        writer.writerows(top10_rows)
    print(f"Done: {summary_path}")
    print(f"Done: {top10_path}")


if __name__ == "__main__":
    main()
