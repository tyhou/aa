# Crossover counts and telomere distance analysis
# Inputs: /data/05_Genetic_Map
# Outputs: /results/05_Genetic_Map

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)

# Paths
input_dir  <- "/data/05_Genetic_Map"
output_dir <- "/results/05_Genetic_Map"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===========================================================================
# Step 1: Generate input files
# ===========================================================================
events <- read.table(file.path(input_dir, "SOLine_15lines_crossover_events_VAFdiff_gt0.8_1391.tsv"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

cat("Original rows:", nrow(events), "\n")

events_filt <- events %>% filter(Abs_Delta_VAF > 0.8)
cat("Rows after VAF > 0.8 filter:", nrow(events_filt), "\n")

# Site-level stats
pos_long <- events_filt %>%
  select(CHROM, Line, SNP1_POS, SNP2_POS) %>%
  pivot_longer(cols = c(SNP1_POS, SNP2_POS), names_to = "SNP_type", values_to = "POS") %>%
  select(CHROM, POS, Line) %>%
  distinct()

pos_counts <- pos_long %>%
  group_by(CHROM, POS) %>%
  summarise(Crossover = n_distinct(Line), .groups = "drop") %>%
  arrange(CHROM, POS)

write.table(pos_counts, file.path(output_dir, "crossover-15lines-pos-only.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Unique sites:", nrow(pos_counts), "\n")

# Window analysis
win_size <- 5000
step     <- 1000

chrom_sizes <- read.table(file.path(input_dir, "genome.chrom.sizes"), header = FALSE,
                          col.names = c("chrom", "length"), stringsAsFactors = FALSE)
chrom_sizes$chrom_chr <- gsub("^Chr", "", as.character(chrom_sizes$chrom))
pos_counts$CHROM_chr   <- gsub("^Chr", "", as.character(pos_counts$CHROM))

window_counts <- pos_counts %>%
  group_by(CHROM) %>%
  group_modify(~ {
    chr_row <- chrom_sizes[chrom_sizes$chrom_chr == gsub("^Chr", "", as.character(.y$CHROM)), ]
    chr_max <- if (nrow(chr_row) > 0) chr_row$length[1] else max(.x$POS, na.rm = TRUE)
    starts  <- seq(1, chr_max, by = step)
    ends    <- starts + win_size - 1
    counts <- sapply(seq_along(starts), function(i) {
      sum(.x$Crossover[.x$POS >= starts[i] & .x$POS <= ends[i]], na.rm = TRUE)
    })
    data.frame(Start = starts, End = ends, Crossover_Count = counts)
  }) %>%
  ungroup()

window_counts <- window_counts %>%
  mutate(CHROM_chr = gsub("^Chr", "", as.character(CHROM))) %>%
  left_join(chrom_sizes %>% select(chrom_chr, length), by = c("CHROM_chr" = "chrom_chr")) %>%
  rename(max_pos = length) %>%
  mutate(End = ifelse(End > max_pos, max_pos, End)) %>%
  select(-max_pos, -CHROM_chr) %>%
  arrange(CHROM, Start) %>%
  select(CHROM, Start, End, Crossover_Count)

write.csv(window_counts, file.path(output_dir, "crossover_window_analysis.csv"), row.names = FALSE)

# ===========================================================================
# Step 2: Heatmap
# ===========================================================================
data <- read.csv(file.path(output_dir, "crossover_window_analysis.csv"), header = TRUE)
colnames(data) <- c("CHROM", "Start", "End", "Crossover_Count")

data$Window_Center_kb <- (data$Start + data$End) / 2 / 1000
data$count_type <- ifelse(data$Crossover_Count == 0, "zero", "non-zero")
data$CHROM <- factor(paste0("Chr", data$CHROM),
                     levels = paste0("Chr", sort(unique(as.numeric(data$CHROM)))))

label_df <- data %>% group_by(CHROM) %>% summarise(x = 0, y = 7, .groups = "drop")
x_max <- max(data$Window_Center_kb, na.rm = TRUE)

p <- ggplot(data, aes(x = Window_Center_kb, y = Crossover_Count)) +
  geom_point(data = subset(data, count_type == "zero"), color = "gray", alpha = 0.5, size = 2.5) +
  geom_point(data = subset(data, count_type == "non-zero"), aes(color = Crossover_Count, alpha = Crossover_Count), size = 3) +
  geom_text(data = label_df, aes(x = x, y = y, label = CHROM), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 6, fontface = "bold") +
  scale_color_gradient(low = "blue", high = "red", name = "Crossover Counts", limits = c(0, 7), oob = squish) +
  scale_alpha_continuous(range = c(0.5, 0.9), guide = "none") +
  scale_x_continuous(labels = comma, limits = c(0, x_max)) +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 7)) +
  facet_wrap(~ CHROM, ncol = 1) +
  labs(x = "Genomic Position (kb)", y = "Crossover Count") +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
        strip.text = element_blank(), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

ggsave(file.path(output_dir, "crossover_points_final.png"), p, width = 13, height = 10, dpi = 300)
ggsave(file.path(output_dir, "crossover_points_final.pdf"), p, width = 13, height = 10, device = "pdf")

# ===========================================================================
# Step 3 & 4: Normalized telomere distance
# ===========================================================================
your_data <- pos_counts
your_data$CHROM_chr <- gsub("^Chr", "", as.character(your_data$CHROM))
your_data$chrom_length <- chrom_sizes$length[match(your_data$CHROM_chr, chrom_sizes$chrom_chr)]

your_data <- your_data[!is.na(your_data$chrom_length), ]

calculate_telomere_distance <- function(position, chrom_length) {
  pmin(position, chrom_length - position) / (chrom_length / 2)
}

your_data$tel_dist <- mapply(calculate_telomere_distance, your_data$POS, your_data$chrom_length)

cat("Mean normalized telomere distance:", round(mean(your_data$tel_dist), 4), "\n")

p1 <- ggplot(your_data, aes(x = tel_dist)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "#C5CAE9", color = "#5C6BC0", alpha = 0.7) +
  geom_density(color = "#B71C1C", linewidth = 1) +
  geom_vline(xintercept = 0.5, color = "grey20", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = mean(your_data$tel_dist), color = "grey20", linewidth = 1) +
  labs(x = "Normalized distance to telomeres", y = "Density") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

p2 <- ggplot(your_data, aes(x = factor(CHROM), y = tel_dist)) +
  geom_boxplot(fill = "#C5CAE9", color = "#5C6BC0", outlier.shape = NA) +
  geom_hline(yintercept = 0.5, color = "grey30", linewidth = 1, linetype = "dashed") +
  geom_point(position = position_jitter(width = 0.2), color = "grey30", alpha = 0.5) +
  scale_x_discrete(labels = function(x) paste0("Chr", gsub("^Chr", "", x))) +
  labs(x = "Chromosome", y = "Normalized distance to telomeres") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

ggsave(file.path(output_dir, "crossover_histogram.png"), p1, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, "crossover_histogram.pdf"), p1, width = 6, height = 4, device = "pdf")
ggsave(file.path(output_dir, "crossover_barplot.png"), p2, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, "crossover_barplot.pdf"), p2, width = 6, height = 4, device = "pdf")

write.csv(your_data, file.path(output_dir, "crossover_analysis_results.csv"), row.names = FALSE)
cat("Crossover counts analysis completed.\n")