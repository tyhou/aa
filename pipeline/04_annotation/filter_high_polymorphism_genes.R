#!/usr/bin/env Rscript

# High-polymorphism gene definition:
#   1) total CDS variants >= 3
#   2) variants per kb CDS >= 1

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg)) setwd(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
summary_file <- if (length(args) >= 1) args[1] else "snpEff_summary.genes.txt"
mapping_file <- if (length(args) >= 2) args[2] else "../../CP01_to_CPHLJ_ID.txt"
gff_file <- if (length(args) >= 3) args[3] else "../../../../../../07-fst-pi-dip/IId_11730_genome.gff3"
output_file <- if (length(args) >= 4) args[4] else "high_polymorphism_genes.xlsx"

min_variants <- 3
min_density <- 1

if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install the R package 'openxlsx'.")
for (f in c(summary_file, mapping_file, gff_file)) if (!file.exists(f)) stop("Input file not found: ", f)

# Read SnpEff gene summary. The first line is a comment and the second line is the header.
snpeff <- read.delim(summary_file, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
names(snpeff)[1] <- sub("^#", "", names(snpeff)[1])

needed <- c("GeneName", "GeneId", "variants_impact_HIGH", "variants_impact_MODERATE",
            "variants_impact_LOW", "variants_effect_missense_variant",
            "variants_effect_synonymous_variant")
if (!all(needed %in% names(snpeff))) stop("Missing required columns in SnpEff summary: ",
                                         paste(setdiff(needed, names(snpeff)), collapse = ", "))

# Map SnpEff CP01 IDs to the CPHLJ reference annotation.
id_map <- read.delim(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
id_map <- id_map[!duplicated(id_map$CP01_ID), ]

# Read GFF3 and extract CDS lengths and functional descriptions by locus_tag.
gff <- read.delim(gff_file, comment.char = "#", header = FALSE, sep = "\t", quote = "",
                  stringsAsFactors = FALSE,
                  col.names = c("Chr", "Source", "Type", "Start", "End",
                                "Score", "Strand", "Phase", "Attributes"))
get_attr <- function(x, key) {
  z <- regmatches(x, regexec(paste0("(?:^|;)", key, "=([^;]*)"), x, perl = TRUE))
  vapply(z, function(a) if (length(a) > 1) utils::URLdecode(a[2]) else NA_character_, character(1))
}
gff$locus_tag <- get_attr(gff$Attributes, "locus_tag")
gff$product <- get_attr(gff$Attributes, "Product")
gff$product[is.na(gff$product)] <- get_attr(gff$Attributes, "product")[is.na(gff$product)]
gff$note <- get_attr(gff$Attributes, "Note")

collapse_text <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(x, collapse = " | ") else NA_character_
}

# Sum the union of CDS intervals, preventing duplicated transcript intervals from double counting.
union_length <- function(start, end) {
  x <- data.frame(start = as.integer(start), end = as.integer(end))
  x <- unique(x[order(x$start, x$end), ])
  if (!nrow(x)) return(NA_integer_)
  total <- 0L; left <- x$start[1]; right <- x$end[1]
  if (nrow(x) > 1) for (i in 2:nrow(x)) {
    if (x$start[i] <= right + 1L) right <- max(right, x$end[i])
    else { total <- total + right - left + 1L; left <- x$start[i]; right <- x$end[i] }
  }
  total + right - left + 1L
}

cds <- gff[gff$Type == "CDS" & !is.na(gff$locus_tag), ]
cds_length <- tapply(seq_len(nrow(cds)), cds$locus_tag,
                     function(i) union_length(cds$Start[i], cds$End[i]))
annotation <- tapply(c(seq_len(nrow(gff))), gff$locus_tag, function(i) {
  p <- collapse_text(gff$product[i]); n <- collapse_text(gff$note[i])
  if (!is.na(p)) p else n
})

x <- merge(snpeff, id_map, by.x = "GeneId", by.y = "CP01_ID", all.x = TRUE, sort = FALSE)
x$CDS_length_bp <- as.integer(cds_length[x$CPHLJ_locus_tag])
x$Annotation <- unname(annotation[x$CPHLJ_locus_tag])
x$Total_CDS_variants <- x$variants_impact_HIGH + x$variants_impact_MODERATE + x$variants_impact_LOW
x$Variant_density_per_kb_CDS <- x$Total_CDS_variants / (x$CDS_length_bp / 1000)

keep <- !is.na(x$CDS_length_bp) & x$CDS_length_bp > 0 &
        x$Total_CDS_variants >= min_variants &
        x$Variant_density_per_kb_CDS >= min_density
out <- x[keep, ]
out <- out[order(-out$Total_CDS_variants, -out$Variant_density_per_kb_CDS, out$GeneId), ]

result <- data.frame(
  No = seq_len(nrow(out)),
  Chr = out$Chromosome,
  `Gene ID (SnpEff CP01)` = out$GeneId,
  `Gene ID (CPHLJ)` = out$CPHLJ_locus_tag,
  Annotation = out$Annotation,
  `No. of variants (CDS)` = out$Total_CDS_variants,
  `Variant density (per kb CDS)` = round(out$Variant_density_per_kb_CDS, 2),
  `CDS length (bp)` = out$CDS_length_bp,
  High = out$variants_impact_HIGH,
  Moderate = out$variants_impact_MODERATE,
  Low = out$variants_impact_LOW,
  Missense = out$variants_effect_missense_variant,
  Synonymous = out$variants_effect_synonymous_variant,
  check.names = FALSE
)

unmatched_index <- is.na(x$CPHLJ_locus_tag) | is.na(x$CDS_length_bp)
unmatched <- data.frame(
  SnpEff_Gene_ID = x$GeneId[unmatched_index],
  Mapping_status = ifelse(is.na(x$CPHLJ_locus_tag[unmatched_index]),
                          "No CP01-to-CPHLJ mapping", "No CDS length in GFF"),
  check.names = FALSE
)

# Create a publication-ready Excel workbook.
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "High-polymorphism genes", gridLines = FALSE)
openxlsx::addWorksheet(wb, "Methods", gridLines = FALSE)
openxlsx::addWorksheet(wb, "Unmatched genes", gridLines = FALSE)
openxlsx::writeData(wb, "High-polymorphism genes", "High-polymorphism genes based on 1,391 variants", startRow = 2)
openxlsx::writeData(wb, "High-polymorphism genes",
                    paste0("Definition: CDS variants >= ", min_variants,
                           " and variant density >= ", min_density, " per kb CDS"), startRow = 3)
openxlsx::writeDataTable(wb, "High-polymorphism genes", result, startRow = 5,
                         tableStyle = "TableStyleMedium2", withFilter = TRUE)

methods <- data.frame(
  Item = c("Input", "CDS variant count", "CDS length", "Variant density", "Filter", "ID note"),
  Definition = c(
    "SnpEff gene summary generated from the VCF containing 1,391 variants.",
    "variants_impact_HIGH + variants_impact_MODERATE + variants_impact_LOW; MODIFIER variants are excluded.",
    "Union length of all CDS intervals assigned to the mapped CPHLJ locus_tag in the GFF3 file.",
    "Total CDS variants / CDS length (bp) x 1,000.",
    "Total CDS variants >= 3 and variant density >= 1 per kb CDS.",
    "The available local mapping links SnpEff CP01 IDs to CPHLJ IDs; no cgd (IIa-IOWA-II) mapping was found."
  ), check.names = FALSE)
openxlsx::writeDataTable(wb, "Methods", methods, tableStyle = "TableStyleMedium2")
openxlsx::writeDataTable(wb, "Unmatched genes", unmatched, tableStyle = "TableStyleMedium2")

title_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold")
openxlsx::addStyle(wb, "High-polymorphism genes", title_style, rows = 2, cols = 1)
openxlsx::freezePane(wb, "High-polymorphism genes", firstActiveRow = 6)
openxlsx::setColWidths(wb, "High-polymorphism genes", cols = 1:13,
                       widths = c(7, 7, 22, 20, 55, 22, 29, 17, 10, 12, 10, 12, 12))
openxlsx::setColWidths(wb, "Methods", cols = 1:2, widths = c(22, 95))
openxlsx::setColWidths(wb, "Unmatched genes", cols = 1:2, widths = c(24, 30))
openxlsx::addStyle(wb, "High-polymorphism genes",
                   openxlsx::createStyle(numFmt = "0.00"), rows = 6:(nrow(result) + 5), cols = 7,
                   gridExpand = TRUE, stack = TRUE)
openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

message("Genes in SnpEff summary: ", nrow(snpeff))
message("High-polymorphism genes: ", nrow(result))
message("Output: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
