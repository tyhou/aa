#!/usr/bin/env Rscript

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

work_dir <- get_script_dir()
input_file <- file.path(work_dir, "snpeff", "tables", "snpeff_primary_effects.tsv")
output_dir <- file.path(work_dir, "SNP_comparison_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dat <- read.delim(input_file, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, check.names = FALSE)

old <- dat[dat$Dataset == "Old_1065" & dat$Variant_Type == "SNP", ]
new <- dat[dat$Dataset == "New_1391" & dat$Variant_Type == "SNP", ]
variant_key <- function(x) paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ":")
old_key <- unique(variant_key(old)); new_key <- unique(variant_key(new))
shared <- intersect(old_key, new_key)
old_only <- setdiff(old_key, new_key)
new_only <- setdiff(new_key, old_key)

if (!identical(c(length(old_key), length(new_key), length(shared), length(old_only), length(new_only)),
               c(1065L, 1138L, 1017L, 48L, 121L))) {
  stop("Unexpected SNP counts; check the SnpEff primary-effects table.")
}

draw_venn <- function() {
  par(mar = c(1, 1, 4, 1), xpd = NA)
  plot(NA, xlim = c(-2.4, 2.4), ylim = c(-1.6, 1.8), asp = 1,
       axes = FALSE, xlab = "", ylab = "",
       main = "Exact SNP overlap: old 1065 vs new 1391 callset")
  symbols(-0.62, 0, circles = 1.35, inches = FALSE, add = TRUE,
          bg = adjustcolor("#4C78A8", alpha.f = 0.38), fg = "#4C78A8", lwd = 2.2)
  symbols( 0.62, 0, circles = 1.35, inches = FALSE, add = TRUE,
          bg = adjustcolor("#E45756", alpha.f = 0.38), fg = "#E45756", lwd = 2.2)
  text(-1.20, 1.12, "Old SNP", cex = 1.25, font = 2)
  text( 1.20, 1.12, "New SNP", cex = 1.25, font = 2)
  text(-1.18, 0.05, length(old_only), cex = 1.8, font = 2)
  text( 0.00, 0.05, length(shared), cex = 1.8, font = 2)
  text( 1.18, 0.05, length(new_only), cex = 1.8, font = 2)
  text(-1.18,-0.22, "old-only", cex = 0.95)
  text( 0.00,-0.22, "shared", cex = 0.95)
  text( 1.18,-0.22, "new-only", cex = 0.95)
  mtext("Exact identity requires identical CHROM, POS, REF and ALT", side = 3, line = 0.5, cex = 0.9)
  text(0, -1.45, "Circle areas are schematic; displayed counts are exact.", cex = 0.78, col = "grey40")
}

png(file.path(output_dir, "SNP_exact_overlap_venn.png"), width = 2400, height = 1800,
    res = 300, bg = "white")
draw_venn(); dev.off()
pdf(file.path(output_dir, "SNP_exact_overlap_venn.pdf"), width = 8, height = 6)
draw_venn(); dev.off()

summary <- data.frame(
  Metric = c("Old_SNP", "New_SNP", "Shared_exact_SNP", "Old_only_SNP", "New_only_SNP"),
  Count = c(length(old_key), length(new_key), length(shared), length(old_only), length(new_only))
)
write.table(summary, file.path(output_dir, "SNP_exact_overlap_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Saved Venn diagram and counts in: ", output_dir)
