# SOLine allele-frequency maps with Heterogeneity.
# Inputs: /data/04_Allele_Frequency/SOLine
# Outputs: /results/04_Allele_Frequency/SOLine

library(ggplot2)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to read Het-SOLine-input.xlsx")
}

# Paths
input_dir  <- "/data/04_Allele_Frequency/SOLine"
output_dir <- "/results/04_Allele_Frequency/SOLine"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Shared helper functions live in this directory; run.sh exports the path
code_dir <- Sys.getenv("ALLELEFREQ_CODE_DIR", unset = "/code/04_Allele_Frequency")

# Source shared functions
source(file.path(code_dir, "allelefreq_functions_1391.R"))

# Read VCF and mask
vcf <- read_haploid_vcf(file.path(input_dir, "SOLine_hap_overlap_1391.vcf"))
mask <- informative_mask(vcf, "11730", "22971")
cat("SOLine VCF records:", nrow(vcf$fixed),
    "\nInformative sites:", sum(mask, na.rm = TRUE), "\n")

freq <- make_frequency_table(vcf, mask, "SOLine_Het_lastDPI", out_dir = output_dir)

# Read metadata and retain max DPI per Line
het_input <- as.data.frame(
  readxl::read_excel(file.path(input_dir, "Het-SOLine-input.xlsx")),
  stringsAsFactors = FALSE
)
required_columns <- c("ID", "Line", "DPI", "Het")
if (!all(required_columns %in% names(het_input))) {
  stop("Het-SOLine-input.xlsx must contain columns: ",
       paste(required_columns, collapse = ", "))
}

het_input$Line_num <- as.integer(sub("^Line\\s*", "", het_input$Line))
het_input$DPI <- as.numeric(het_input$DPI)
het_input$Het <- as.numeric(het_input$Het)

max_dpi <- ave(het_input$DPI, het_input$Line_num, FUN = max)
sample_meta <- het_input[het_input$DPI == max_dpi, , drop = FALSE]
sample_meta <- sample_meta[order(sample_meta$Line_num), , drop = FALSE]

if (anyDuplicated(sample_meta$Line_num)) {
  stop("More than one sample is tied for the maximum DPI in at least one Line")
}
if (any(!is.finite(sample_meta$Het))) {
  stop("Missing or non-numeric Het value in the selected last-DPI samples")
}

sample_meta <- data.frame(
  WGA = sample_meta$ID,
  Line = sample_meta$Line_num,
  DPI = sample_meta$DPI,
  Sample_info = sprintf("Line%d_D%d", sample_meta$Line_num, as.integer(sample_meta$DPI)),
  Plot_label = sprintf("Line%d_D%d (Heterogeneity= %.1f)",
                       sample_meta$Line_num, as.integer(sample_meta$DPI), sample_meta$Het),
  Het = sample_meta$Het,
  stringsAsFactors = FALSE
)

missing_samples <- setdiff(sample_meta$WGA, vcf$samples)
if (length(missing_samples)) {
  stop("Selected last-DPI samples missing from SOLine_hap_overlap_1391.vcf: ",
       paste(missing_samples, collapse = ", "))
}

cat("Selected last-DPI samples by Line:\n")
print(sample_meta[, c("WGA", "Line", "DPI", "Het", "Plot_label")], row.names = FALSE)

# Define labels
labels <- setNames(vcf$samples, vcf$samples)
labels[sample_meta$WGA] <- sample_meta$Plot_label
labels[c("22971", "11730")] <- c("GD", "HLJ")

# Define groups
soline_groups <- list(
  SOLine_15lines_Het_lastDPI = sample_meta$WGA,
  SOLine_GD_HLJ_blank_15lines_Het_lastDPI = c(
    "22971", "11730", NA_character_, sample_meta$WGA
  )
)

# Plot maps
plot_frequency_maps(
  freq, vcf$contigs, "SOLine", "SOLine_Het_lastDPI", soline_groups, labels,
  out_dir = output_dir
)

# Calculate crossovers
crossover <- calculate_crossovers(
  vcf, mask, sample_meta, "SOLine_15lines_Het_lastDPI", threshold = 0.8,
  out_dir = output_dir
)

print(crossover$summary)
cat("Outputs written to:", output_dir, "\n")
cat("SOLine Het-labelled last-DPI analysis completed.\n")