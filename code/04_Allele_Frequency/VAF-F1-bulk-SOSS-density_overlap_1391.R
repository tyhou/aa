# VAF density analysis
# Inputs: /data/04_Allele_Frequency/VAF_Density
# Outputs: /results/04_Allele_Frequency/VAF_Density

library(ggplot2)
library(dplyr)

# Paths
input_dir  <- "/data/04_Allele_Frequency/VAF_Density"
output_dir <- "/results/04_Allele_Frequency/VAF_Density"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define groups
parents <- c("11730", "22971")
f1_pool <- c("W347")
oocyst <- c("W585", "W586", "W587", "W588", "W589", "W590", 
            "W591", "W592", "W594", "W595", "W596")
sporozoite <- c("W680", "W681", "W682",
                "W788", "W789", "W790", "W791",
                "W792", "W793", "W794", "W795")
samples_of_interest <- c(parents, f1_pool, oocyst, sporozoite)

sample_colors <- c(
  "11730" = "#DA2C2C", "22971" = "#2C9B2C", "W347"  = "#FDD835",
  "W585"  = "#FDD835", "W586"  = "#FDD835", "W587"  = "#FDD835",
  "W588"  = "#FDD835", "W589"  = "#FDD835", "W590"  = "#FDD835",
  "W591"  = "#FDD835", "W592"  = "#FDD835", "W594"  = "#FDD835",
  "W595"  = "#FDD835", "W596"  = "#FDD835",
  "W680"  = "#1E88E5", "W681"  = "#1E88E5", "W682"  = "#1E88E5",
  "W788"  = "#1E88E5", "W789"  = "#1E88E5", "W790"  = "#1E88E5",
  "W791"  = "#1E88E5", "W792"  = "#1E88E5", "W793"  = "#1E88E5",
  "W794"  = "#1E88E5", "W795"  = "#1E88E5"
)

# Read VCF
vcf_file <- file.path(input_dir, "F1-bulk-SOSS_dip_overlap_1391.vcf")
vcf_lines <- readLines(vcf_file)
header_line <- grep("^#CHROM", vcf_lines)
if (length(header_line) == 0) stop("Cannot find #CHROM header line")

header <- strsplit(vcf_lines[header_line], "\t")[[1]]
sample_cols <- header[10:length(header)]
data_lines <- vcf_lines[(header_line + 1):length(vcf_lines)]
data_lines <- data_lines[data_lines != ""]

# Parse GT and AD
vaf_list <- list()
samples_to_process <- intersect(sample_cols, samples_of_interest)
missing_samples <- setdiff(samples_of_interest, sample_cols)
if (length(missing_samples) > 0) warning("Missing samples: ", paste(missing_samples, collapse = ", "))

for (sample in samples_to_process) {
  sample_idx <- match(sample, header)
  vaf_values <- numeric(0)
  
  for (line in data_lines) {
    cols <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (sample_idx > length(cols) || length(cols) < 10) next
    
    format_keys <- strsplit(cols[9], ":", fixed = TRUE)[[1]]
    sample_values <- strsplit(cols[sample_idx], ":", fixed = TRUE)[[1]]
    gt_idx <- match("GT", format_keys)
    ad_idx <- match("AD", format_keys)
    if (is.na(gt_idx) || is.na(ad_idx) || gt_idx > length(sample_values) || ad_idx > length(sample_values)) next
    
    if (sample_values[gt_idx] %in% c(".", "./.", ".|.")) next
    
    ad_values <- suppressWarnings(as.numeric(strsplit(sample_values[ad_idx], ",", fixed = TRUE)[[1]]))
    if (length(ad_values) >= 2 && !anyNA(ad_values) && sum(ad_values) > 0) {
      vaf_values <- c(vaf_values, sum(ad_values[-1]) / sum(ad_values))
    }
  }
  
  if (length(vaf_values) > 0) {
    vaf_list[[sample]] <- data.frame(Sample = sample, VAF = vaf_values)
  }
}

vaf_results <- do.call(rbind, vaf_list)
total_snps <- length(data_lines)

# Add labels and colors
vaf_results$Type <- "Other"
vaf_results$Type[vaf_results$Sample %in% parents] <- "Parent"
vaf_results$Type[vaf_results$Sample %in% f1_pool] <- "F1_Pool"
vaf_results$Type[vaf_results$Sample %in% oocyst] <- "Single_Oocyst"
vaf_results$Type[vaf_results$Sample %in% sporozoite] <- "Single_Sporozoite"
vaf_results$Color <- sample_colors[vaf_results$Sample]

# Calculate VAF statistics
vaf_stats <- vaf_results %>%
  group_by(Sample, Type) %>%
  summarise(
    n_sites = n(),
    pct_of_total = n() / total_snps * 100,
    mean_vaf = mean(VAF, na.rm = TRUE),
    median_vaf = median(VAF, na.rm = TRUE),
    pct_0_10 = sum(VAF >= 0 & VAF < 0.1) / n() * 100,
    pct_10_20 = sum(VAF >= 0.1 & VAF < 0.2) / n() * 100,
    pct_20_30 = sum(VAF >= 0.2 & VAF < 0.3) / n() * 100,
    pct_30_40 = sum(VAF >= 0.3 & VAF < 0.4) / n() * 100,
    pct_40_50 = sum(VAF >= 0.4 & VAF < 0.5) / n() * 100,
    pct_50_60 = sum(VAF >= 0.5 & VAF < 0.6) / n() * 100,
    pct_60_70 = sum(VAF >= 0.6 & VAF < 0.7) / n() * 100,
    pct_70_80 = sum(VAF >= 0.7 & VAF < 0.8) / n() * 100,
    pct_80_90 = sum(VAF >= 0.8 & VAF < 0.9) / n() * 100,
    pct_90_100 = sum(VAF >= 0.9 & VAF <= 1) / n() * 100,
    .groups = "drop"
  )

write.table(vaf_stats, file = file.path(output_dir, "VAF_Statistics_overlap_1391.xls"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plot VAF distribution
bin_width <- 0.02
bins <- seq(0, 1, by = bin_width)
bins[length(bins)] <- 1 + .Machine$double.eps

bin_data <- vaf_results %>%
  mutate(bin = cut(VAF, breaks = bins, include.lowest = TRUE, right = FALSE)) %>%
  filter(!is.na(bin)) %>%
  group_by(Sample, Type, Color, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(pct_of_sites = count / sum(count) * 100) %>%
  ungroup()

bin_midpoints <- seq(bin_width/2, 1 - bin_width/2, by = bin_width)
bin_data$bin_mid <- bin_midpoints[as.numeric(bin_data$bin)]

sample_order <- c("11730", "22971", "W347",
                  "W585", "W586", "W587", "W588", "W589", "W590",
                  "W591", "W592", "W594", "W595", "W596",
                  "W680", "W681", "W682",
                  "W788", "W789", "W790", "W791",
                  "W792", "W793", "W794", "W795")
bin_data$Sample <- factor(bin_data$Sample, levels = sample_order)

p <- ggplot(bin_data, aes(x = bin_mid, y = pct_of_sites, fill = Color)) +
  geom_col(width = bin_width * 0.9) +
  scale_fill_identity() +
  facet_wrap(~Sample, ncol = 5) +
  labs(
    title = "VAF Distribution by Sample",
    subtitle = "Parent (red/green) | F1 Pool (yellow) | Single Oocyst (yellow) | Single Sporozoite (blue)",
    x = "Variant Allele Frequency (VAF)",
    y = "Percentage of Valid VAF Sites (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 0.5
  )

ggsave(file.path(output_dir, "VAF_Distribution_All_Samples_overlap_1391.pdf"), plot = p, width = 18, height = 15)
ggsave(file.path(output_dir, "VAF_Distribution_All_Samples_overlap_1391.png"), plot = p, width = 18, height = 15, dpi = 300)

# Output summary
cat("VAF Analysis Complete.\nTotal SNPs:", total_snps, 
    "\nTotal VAF values:", nrow(vaf_results), 
    "\nSamples analyzed:", length(unique(vaf_results$Sample)), "\n")