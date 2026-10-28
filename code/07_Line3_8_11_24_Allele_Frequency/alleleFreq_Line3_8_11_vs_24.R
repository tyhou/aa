#!/usr/bin/env Rscript

# Lines 3/8/11 (High) vs Line24 (Low)
# Candidate threshold: BH-FDR < 0.05 and VAF_High - VAF_Low < -0.5

# Paths
INPUT_DIR  <- "/data/07_Line3_8_11_24_Allele_Frequency"
OUTPUT_DIR <- "/results/07_Line3_8_11_24_Allele_Frequency"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

vaf_file <- file.path(INPUT_DIR, "SOLine_alt_allele_frequency_1391.tsv")
group_file <- file.path(INPUT_DIR, "selected_samples.tsv")
gff_file <- file.path(INPUT_DIR, "IId_11730_genome.gff3")

fdr_cutoff <- 0.05
effect_cutoff <- -0.5
min_nonmissing_per_group <- 2L

for (f in c(vaf_file, group_file, gff_file)) if (!file.exists(f)) stop("Missing input: ", f)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx")

# Read metadata
meta <- read.delim(group_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("vcf_id", "Group", "Analysis_group") %in% names(meta)))
  stop("selected_samples.tsv needs vcf_id, Group, and Analysis_group columns")
if (!setequal(unique(meta$Analysis_group), c("High", "Low"))) stop("Analysis_group must contain High and Low")

# Read VAF data
dat <- read.delim(vaf_file, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("NA", "NaN"), quote = "")
variant_cols <- c("CHROM", "POS", "ID", "REF", "ALT")
if (!all(c(variant_cols, meta$vcf_id) %in% names(dat))) stop("VAF table lacks required columns or samples")
for (s in meta$vcf_id) dat[[s]] <- suppressWarnings(as.numeric(dat[[s]]))

high_ids <- meta$vcf_id[meta$Analysis_group == "High"]
low_ids <- meta$vcf_id[meta$Analysis_group == "Low"]
high_mat <- as.matrix(dat[, high_ids, drop = FALSE])
low_mat <- as.matrix(dat[, low_ids, drop = FALSE])

# Calculate weighted means
high_line_ids <- split(meta$vcf_id[meta$Analysis_group == "High"], meta$Group[meta$Analysis_group == "High"])
high_line_means <- sapply(high_line_ids, function(ids) rowMeans(dat[, ids, drop = FALSE], na.rm = TRUE))
mean_high <- rowMeans(high_line_means, na.rm = TRUE)
mean_low <- rowMeans(low_mat, na.rm = TRUE)
mean_high[!is.finite(mean_high)] <- NA_real_
mean_low[!is.finite(mean_low)] <- NA_real_

# Wilcoxon test
wilcox_p <- vapply(seq_len(nrow(dat)), function(i) {
  x <- high_mat[i, ]; y <- low_mat[i, ]
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < min_nonmissing_per_group || length(y) < min_nonmissing_per_group) return(NA_real_)
  if (length(unique(c(x, y))) < 2L) return(1)
  suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
}, numeric(1))

result <- data.frame(dat[, variant_cols], n_High = rowSums(!is.na(high_mat)),
  n_Low = rowSums(!is.na(low_mat)), mean_VAF_High = mean_high, mean_VAF_Low = mean_low,
  VAF_difference_High_minus_Low = mean_high - mean_low, wilcox_p = wilcox_p,
  FDR = p.adjust(wilcox_p, method = "BH"), check.names = FALSE)
result$candidate <- !is.na(result$FDR) & result$FDR < fdr_cutoff &
  !is.na(result$VAF_difference_High_minus_Low) & result$VAF_difference_High_minus_Low < effect_cutoff

# Parse GFF3 and annotate
gff <- read.delim(gff_file, comment.char = "#", header = FALSE, sep = "\t", quote = "",
  stringsAsFactors = FALSE, col.names = c("seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"))
attr_value <- function(x, key) {
  z <- regmatches(x, regexec(paste0("(?:^|;)", key, "=([^;]*)"), x, perl = TRUE))
  vapply(z, function(a) if (length(a) > 1L) utils::URLdecode(a[2]) else NA_character_, character(1))
}
gff$locus_tag <- attr_value(gff$attributes, "locus_tag")
gff$product <- attr_value(gff$attributes, "Product")
gff$product[is.na(gff$product)] <- attr_value(gff$attributes, "product")[is.na(gff$product)]
gff$note <- attr_value(gff$attributes, "Note")
genes <- gff[gff$type == "gene", c("seqid", "start", "end", "strand", "locus_tag")]
genes$gene_id <- sub("^gene-", "", attr_value(gff$attributes[gff$type == "gene"], "ID"))
collapse_text <- function(x) { x <- unique(x[!is.na(x) & nzchar(x)]); if (length(x)) paste(x, collapse = " | ") else NA_character_ }
product_map <- tapply(gff$product, gff$locus_tag, collapse_text)
note_map <- tapply(gff$note, gff$locus_tag, collapse_text)
genes$product <- unname(product_map[genes$locus_tag]); genes$note <- unname(note_map[genes$locus_tag])
annotate_site <- function(chr, pos) {
  z <- genes[as.character(genes$seqid) == as.character(chr), , drop = FALSE]
  overlap <- z$start <= pos & z$end >= pos
  if (any(overlap)) { z <- z[overlap, , drop = FALSE]; relation <- "overlap"; distance <- 0L }
  else { d <- pmin(abs(pos - z$start), abs(pos - z$end)); z <- z[d == min(d), , drop = FALSE]; relation <- "nearest"; distance <- min(d) }
  data.frame(gene_relation = relation, distance_to_gene = distance,
    gene_id = collapse_text(z$gene_id), locus_tag = collapse_text(z$locus_tag),
    gene_start = paste(z$start, collapse = " | "), gene_end = paste(z$end, collapse = " | "),
    strand = collapse_text(z$strand), product = collapse_text(z$product), note = collapse_text(z$note))
}
candidates <- result[result$candidate, , drop = FALSE]
candidates <- cbind(candidates, do.call(rbind, Map(annotate_site, candidates$CHROM, candidates$POS)))
candidates <- candidates[order(candidates$VAF_difference_High_minus_Low, candidates$FDR), ]
write.table(candidates, file.path(OUTPUT_DIR, "candidate_VAF_sites_annotated.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

# Plot volcano
library(ggplot2)
result$Category <- ifelse(result$candidate, "Significant", "Others")
p1 <- ggplot(result, aes(VAF_difference_High_minus_Low, -log10(pmax(FDR, 1e-300)), colour = Category)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", colour = "#E60012") +
  geom_vline(xintercept = effect_cutoff, linetype = "dashed", colour = "#E60012") +
  scale_colour_manual(values = c(Significant = "#E60012", Others = "#888888"),
    labels = c(Significant = "Significant (FDR < 0.05 & Effect size < -0.5)", Others = "Others")) +
  labs(x = "Effect size (VAF_High - VAF_Low)", y = expression(-log[10]~FDR), colour = NULL) +
  theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), legend.position = "right")

# Plot boxplot for group means
group_mean_vaf <- rbind(
  data.frame(group = "High", VAF = candidates$mean_VAF_High),
  data.frame(group = "Low", VAF = candidates$mean_VAF_Low)
)
group_mean_vaf$group <- factor(group_mean_vaf$group, c("High", "Low"))
set.seed(20260908)
p2 <- ggplot(group_mean_vaf, aes(group, VAF, colour = group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 1.5) +
  scale_colour_manual(values = c(High = "#E60012", Low = "#377EB8")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Variant Allele Frequency (VAF)") +
  theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), legend.position = "none")

# Save combined plots
plot_png <- file.path(OUTPUT_DIR, "candidate_VAF_plot.png")
png(plot_png, width = 1800, height = 760, res = 160)
grid::grid.newpage(); print(p1, vp = grid::viewport(x = 0.35, width = 0.68)); print(p2, vp = grid::viewport(x = 0.84, width = 0.30)); dev.off()

pdf(file.path(OUTPUT_DIR, "candidate_VAF_plot.pdf"), width = 14, height = 3.5, useDingbats = FALSE)
grid::grid.newpage(); print(p1, vp = grid::viewport(x = 0.35, width = 0.68)); print(p2, vp = grid::viewport(x = 0.84, width = 0.30)); dev.off()

# Save Excel report
candidate_genes <- unique(candidates[, c("CHROM", "POS", "REF", "ALT", "mean_VAF_High",
  "mean_VAF_Low", "VAF_difference_High_minus_Low", "FDR", "gene_relation",
  "distance_to_gene", "gene_id", "locus_tag", "gene_start", "gene_end", "strand", "product", "note")])
wb <- openxlsx::createWorkbook()
for (s in c("Summary", "Candidate_sites", "Candidate_genes", "Samples")) openxlsx::addWorksheet(wb, s, gridLines = FALSE)
summary_data <- data.frame(Metric = c("Comparison", "Candidate threshold", "Variants tested", "Candidate sites"),
  Value = c("Lines 3/8/11 High vs Line24 Low", "FDR < 0.05 and VAF_High - VAF_Low < -0.5", nrow(result), nrow(candidates)))
openxlsx::writeDataTable(wb, "Summary", summary_data, tableStyle = "TableStyleMedium2")
openxlsx::insertImage(wb, "Summary", plot_png, startRow = 8, startCol = 1, width = 12, height = 5.1)
openxlsx::writeDataTable(wb, "Candidate_sites", candidates, tableStyle = "TableStyleMedium2")
openxlsx::writeDataTable(wb, "Candidate_genes", candidate_genes, tableStyle = "TableStyleMedium2")
openxlsx::writeDataTable(wb, "Samples", meta, tableStyle = "TableStyleMedium2")
for (s in c("Candidate_sites", "Candidate_genes", "Samples")) openxlsx::freezePane(wb, s, firstRow = TRUE)
openxlsx::setColWidths(wb, "Summary", 1:2, c(24, 55))
openxlsx::saveWorkbook(wb, file.path(OUTPUT_DIR, "candidate_VAF_sites_and_genes.xlsx"), overwrite = TRUE)

# Save summary text
writeLines(c(paste0("High: Lines 3/8/11 (n=", length(high_ids), "); Low: Line24 (n=", length(low_ids), ")"),
  "Candidate threshold: FDR < 0.05 and VAF_High - VAF_Low < -0.5",
  paste0("Variants tested: ", nrow(result)), paste0("Candidate sites: ", nrow(candidates)),
  "Wilcoxon/FDR is exploratory because samples from the same line are repeated observations."),
  file.path(OUTPUT_DIR, "candidate_VAF_summary.txt"))
message("Candidate sites: ", nrow(candidates))
cat("Mean VAF High (candidates):", mean(candidates$mean_VAF_High, na.rm = TRUE), "\n")
cat("Mean VAF Low (candidates):",  mean(candidates$mean_VAF_Low,  na.rm = TRUE), "\n")