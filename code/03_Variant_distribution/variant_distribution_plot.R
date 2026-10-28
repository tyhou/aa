# Load packages
library(ggplot2)
library(dplyr)
library(patchwork)

# Paths
input_dir  <- "/data/03_Variant_distribution"
output_dir <- "/results/03_Variant_distribution"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Chromosome lengths
chrom_lengths <- data.frame(
  CHROM = c("1", "2", "3", "4", "5", "6", "7", "8"),
  length = c(883282, 990770, 1100651, 1105975, 1090562, 1300158, 1327643, 1332802)
)

# Read positions
snp_positions <- read.table(file.path(input_dir, "intersection_snps_positions.txt"), 
                             header = FALSE, stringsAsFactors = FALSE)
colnames(snp_positions) <- c("CHROM", "POS")
snp_positions$CHROM <- as.character(snp_positions$CHROM)

indel_positions <- read.table(file.path(input_dir, "intersection_indels_positions.txt"), 
                               header = FALSE, stringsAsFactors = FALSE)
colnames(indel_positions) <- c("CHROM", "POS")
indel_positions$CHROM <- as.character(indel_positions$CHROM)

all_positions <- rbind(
  cbind(snp_positions, TYPE = "SNP"),
  cbind(indel_positions, TYPE = "InDel")
)

# Calculate window counts
window_size <- 1000

calculate_window_counts <- function(variant_data, chrom_lengths, window_size) {
  complete_data <- data.frame()
  
  for(i in 1:nrow(chrom_lengths)) {
    chrom <- chrom_lengths$CHROM[i]
    total_bins <- ceiling(chrom_lengths$length[i] / window_size)
    
    all_bins <- data.frame(
      CHROM = chrom,
      BIN_1Kb = 0:(total_bins-1)
    )
    
    chrom_variants <- variant_data[variant_data$CHROM == chrom, ]
    if(nrow(chrom_variants) > 0) {
      chrom_variants$BIN_1Kb <- floor(chrom_variants$POS / window_size)
      variant_counts <- chrom_variants %>% 
        group_by(BIN_1Kb) %>% 
        summarise(VARIANTS = n(), .groups = 'drop')
    } else {
      variant_counts <- data.frame(BIN_1Kb = integer(0), VARIANTS = integer(0))
    }
    
    chrom_complete <- merge(all_bins, variant_counts, by = "BIN_1Kb", all.x = TRUE)
    chrom_complete$VARIANTS[is.na(chrom_complete$VARIANTS)] <- 0
    chrom_complete$CHROM <- chrom
    
    complete_data <- rbind(complete_data, chrom_complete)
  }
  
  complete_data <- complete_data[order(complete_data$CHROM, complete_data$BIN_1Kb), ]
  return(complete_data)
}

# Calculate for datasets
data_snp <- calculate_window_counts(snp_positions, chrom_lengths, window_size)
data_indel <- calculate_window_counts(indel_positions, chrom_lengths, window_size)
data_all <- calculate_window_counts(all_positions, chrom_lengths, window_size)

# Add labels
data_snp$DATASET <- paste0("SNP (n=", nrow(snp_positions), ")")
data_indel$DATASET <- paste0("InDel (n=", nrow(indel_positions), ")")
data_all$DATASET <- paste0("All Variants (n=", nrow(all_positions), ")")

# Cumulative positions
add_cumulative_positions <- function(data, chrom_lengths, window_size) {
  data$POS_Mb <- data$BIN_1Kb * window_size / 1000000
  data$cumulative_pos <- 0
  
  cumulative_pos <- 0
  for(chrom in unique(chrom_lengths$CHROM)) {
    chrom_data <- data[data$CHROM == chrom, ]
    if(nrow(chrom_data) > 0) {
      max_pos <- max(chrom_data$POS_Mb)
      data$cumulative_pos[data$CHROM == chrom] <- 
        cumulative_pos + data$POS_Mb[data$CHROM == chrom]
      cumulative_pos <- cumulative_pos + max_pos + 0.1
    }
  }
  
  chromosome_breaks <- c()
  chromosome_labels <- c()
  cumulative_pos <- 0
  for(chrom in unique(chrom_lengths$CHROM)) {
    max_pos <- max(data$POS_Mb[data$CHROM == chrom])
    chromosome_breaks <- c(chromosome_breaks, cumulative_pos + max_pos/2)
    chromosome_labels <- c(chromosome_labels, paste0("Chr", chrom))
    cumulative_pos <- cumulative_pos + max_pos + 0.1
  }
  
  return(list(data = data, breaks = chromosome_breaks, labels = chromosome_labels))
}

result_snp <- add_cumulative_positions(data_snp, chrom_lengths, window_size)
result_indel <- add_cumulative_positions(data_indel, chrom_lengths, window_size)
result_all <- add_cumulative_positions(data_all, chrom_lengths, window_size)

data_snp <- result_snp$data
data_indel <- result_indel$data
data_all <- result_all$data

common_breaks <- result_all$breaks
common_labels <- result_all$labels

# Colors
data_snp$chrom_color <- ifelse(as.numeric(data_snp$CHROM) %% 2 == 1, "gray50", "black")
data_indel$chrom_color <- ifelse(as.numeric(data_indel$CHROM) %% 2 == 1, "gray50", "black")
data_all$chrom_color <- ifelse(as.numeric(data_all$CHROM) %% 2 == 1, "gray50", "black")

# Plots
p1 <- ggplot(data_snp, aes(x = cumulative_pos, y = VARIANTS)) +
  geom_col(aes(color = chrom_color), fill = "#1f78b4", linewidth = 0.3) +
  scale_color_identity() +
  scale_x_continuous(breaks = common_breaks, labels = common_labels) +
  labs(title = "", x = "", y = "Count per 1kb") +
  ylim(0, max(data_all$VARIANTS) * 1.1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

p2 <- ggplot(data_indel, aes(x = cumulative_pos, y = VARIANTS)) +
  geom_col(aes(color = chrom_color), fill = "#e31a1c", linewidth = 0.3) +
  scale_color_identity() +
  scale_x_continuous(breaks = common_breaks, labels = common_labels) +
  labs(title = "", x = "", y = "Count per 1kb") +
  ylim(0, max(data_all$VARIANTS) * 1.1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

p3 <- ggplot(data_all, aes(x = cumulative_pos, y = VARIANTS)) +
  geom_col(aes(color = chrom_color), fill = "#33a02c", linewidth = 0.3) +
  scale_color_identity() +
  scale_x_continuous(breaks = common_breaks, labels = common_labels) +
  labs(title = "Variant Density Distribution (1kb Windows)", x = "Chromosome", y = "Count per 1kb") +
  ylim(0, max(data_all$VARIANTS) * 1.1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Combine
combined_plot <- (p1 + labs(title = "SNP")) / 
                 (p2 + labs(title = "InDel")) / 
                 (p3 + labs(title = "All Variants")) +
  plot_annotation(
    title = "Variant Distribution Across 8 Chromosomes",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  ) &
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# Save plots
ggsave(file.path(output_dir, "variant_distribution_combined.png"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "variant_distribution_combined.pdf"), combined_plot, width = 12, height = 10, device = "pdf")

ggsave(file.path(output_dir, "SNP_distribution.png"), p1, width = 12, height = 3.5, dpi = 300)
ggsave(file.path(output_dir, "InDel_distribution.png"), p2, width = 12, height = 3.5, dpi = 300)
ggsave(file.path(output_dir, "All_variants_distribution.png"), p3, width = 12, height = 3.5, dpi = 300)

# Print stats
cat("\n========================================\n")
cat("Variant Statistics Summary\n")
cat("========================================\n\n")

cat("SNP total:", nrow(snp_positions), "\n")
cat("InDel total:", nrow(indel_positions), "\n")
cat("All variants total:", nrow(all_positions), "\n")

cat("\n--- SNP distribution by chromosome ---\n")
print(table(snp_positions$CHROM))

cat("\n--- InDel distribution by chromosome ---\n")
print(table(indel_positions$CHROM))

cat("\n--- All variants distribution by chromosome ---\n")
print(table(all_positions$CHROM))

cat("\n========================================\n")
cat("Plots saved:\n")
cat("  - variant_distribution_combined.png/pdf\n")
cat("  - SNP_distribution.png\n")
cat("  - InDel_distribution.png\n")
cat("  - All_variants_distribution.png\n")
cat("========================================\n")

print(combined_plot)