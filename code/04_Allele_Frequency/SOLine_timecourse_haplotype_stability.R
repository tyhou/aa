# SOLine time-course haplotype stability.
# Inputs: /data/04_Allele_Frequency/SOLine
# Outputs: /results/04_Allele_Frequency/SOLine

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# Paths
input_dir  <- "/data/04_Allele_Frequency/SOLine"
output_dir <- "/results/04_Allele_Frequency/SOLine"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read VAF table
vaf_file <- file.path(input_dir, "SOLine_alt_allele_frequency_1391.tsv")
vaf <- read.delim(vaf_file, stringsAsFactors = FALSE, check.names = FALSE)

# Define sample pairs
line_pairs <- data.frame(
  Line  = c("Line1","Line2","Line3","Line4","Line6","Line7","Line8",
            "Line10","Line11","Line16","Line19","Line24","Line30","Line31","Line32"),
  Early = c("W683","W685","W373","W688","W469","W376","W447",
            "W377","W378","W379","W380","W441","W388","W383","W384"),
  Late  = c("W455","W454","W468","W689","W375","W456","W453",
            "W457","W470","W458","W459","W446","W464","W465","W466"),
  stringsAsFactors = FALSE
)

# Order lines by number, then reverse
line_pairs$Line_num <- as.numeric(gsub("Line", "", line_pairs$Line))
line_pairs <- line_pairs[order(line_pairs$Line_num), ]
line_levels <- rev(line_pairs$Line)

all_samples <- colnames(vaf)
missing <- setdiff(c(line_pairs$Early, line_pairs$Late), all_samples)
if (length(missing) > 0) {
  warning("Missing samples in VAF table: ", paste(missing, collapse = ", "))
}

# Build long table per line
line_results <- list()

for (i in seq_len(nrow(line_pairs))) {
  ln  <- line_pairs$Line[i]
  er  <- line_pairs$Early[i]
  lt  <- line_pairs$Late[i]
  if (!(er %in% all_samples) || !(lt %in% all_samples)) next
  
  df <- vaf %>%
    select(CHROM, POS, all_of(er), all_of(lt)) %>%
    rename(Early_VAF = all_of(er), Late_VAF = all_of(lt)) %>%
    mutate(
      CHROM      = factor(CHROM, levels = sort(unique(as.numeric(as.character(CHROM))))),
      Early_VAF  = as.numeric(Early_VAF),
      Late_VAF   = as.numeric(Late_VAF)
    ) %>%
    filter(!is.na(Early_VAF) & !is.na(Late_VAF)) %>%
    mutate(
      Line          = ln,
      Delta_VAF     = Late_VAF - Early_VAF,
      Abs_Delta_VAF = abs(Delta_VAF),
      Informative   = (Early_VAF > 0.15 & Early_VAF < 0.85) |
        (Late_VAF  > 0.15 & Late_VAF  < 0.85)
    )
  
  if (nrow(df) > 0) line_results[[ln]] <- df
}

all_df <- bind_rows(line_results)

# Stability stats per line
summary_stats <- all_df %>%
  group_by(Line) %>%
  summarise(
    n_snps          = n(),
    n_informative   = sum(Informative, na.rm = TRUE),
    pearson_r       = cor(Early_VAF, Late_VAF, use = "complete.obs"),
    spearman_rho    = cor(Early_VAF, Late_VAF, method = "spearman", use = "complete.obs"),
    median_abs_delta= median(Abs_Delta_VAF, na.rm = TRUE),
    p95_abs_delta   = quantile(Abs_Delta_VAF, 0.95, na.rm = TRUE),
    frac_stable     = mean(Abs_Delta_VAF < 0.2, na.rm = TRUE),
    frac_shift_gt02 = mean(Abs_Delta_VAF > 0.2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Line)

print(summary_stats)

# Histogram of |ΔVAF|
p_hist <- ggplot(all_df, aes(x = Abs_Delta_VAF)) +
  geom_histogram(bins = 50, fill = "tomato", color = "white") +
  geom_vline(xintercept = 0.2, linetype = 2, color = "blue") +
  geom_vline(xintercept = 0.5, linetype = 2, color = "darkblue") +
  facet_wrap(~ Line, ncol = 5, scales = "free_y") +
  labs(
    x = "|VAF(late) - VAF(early)|",
    y = "SNP count",
    title = "Allele-frequency change between early and late time points",
    subtitle = "Blue dashed lines = 0.2 and 0.5 thresholds"
  ) +
  theme_bw(base_size = 9)

# Pearson correlation per line
p_corr <- ggplot(summary_stats, aes(x = Line, y = pearson_r)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(pearson_r, 3)), vjust = -0.3, size = 2.8) +
  ylim(0, 1.05) +
  labs(x = NULL, y = "Pearson r (early vs late VAF)",
       title = "VAF correlation between earliest and latest time points") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Fraction of stable SNPs
p_stable <- ggplot(summary_stats, aes(x = Line, y = frac_stable)) +
  geom_col(fill = "forestgreen") +
  geom_text(aes(label = scales::percent(frac_stable, accuracy = 1)),
            vjust = -0.3, size = 2.8) +
  ylim(0, 1.05) +
  labs(x = NULL, y = "Fraction of SNPs with |ΔVAF| < 0.2",
       title = "Fraction of stable SNPs between time points") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate crossover breakpoints
compute_crossovers <- function(df, vaf_col, threshold = 0.8) {
  df <- df %>% arrange(CHROM, POS)
  out <- list()
  for (chr in levels(df$CHROM)) {
    sub <- df[df$CHROM == chr, ]
    if (nrow(sub) < 2) next
    v <- sub[[vaf_col]]
    d <- diff(v)
    idx <- which(abs(d) > threshold)
    if (length(idx) > 0) {
      out[[chr]] <- data.frame(
        CHROM = chr,
        POS1  = sub$POS[idx],
        POS2  = sub$POS[idx + 1],
        Delta = d[idx]
      )
    }
  }
  bind_rows(out)
}

cross_list <- list()
for (ln in names(line_results)) {
  df <- line_results[[ln]]
  for (tp in c("Early", "Late")) {
    col <- if (tp == "Early") "Early_VAF" else "Late_VAF"
    cr <- compute_crossovers(df, col)
    if (nrow(cr) > 0) {
      cr$Line <- ln
      cr$Timepoint <- tp
      cross_list[[paste(ln, tp, sep = "_")]] <- cr
    }
  }
}
cross_df <- bind_rows(cross_list)

cross_df$Line <- factor(cross_df$Line, levels = line_levels)

cross_df$CHROM_label <- factor(paste0("Chr", cross_df$CHROM),
                               levels = paste0("Chr", sort(unique(as.numeric(as.character(cross_df$CHROM))))))

if (nrow(cross_df) > 0) {
  
  # Identify shared breakpoints
  early_df <- cross_df %>% filter(Timepoint == "Early") %>%
    mutate(key = paste(Line, CHROM, POS1, sep = "_")) %>%
    select(Line, CHROM, CHROM_label, POS1, key)
  
  late_df <- cross_df %>% filter(Timepoint == "Late") %>%
    mutate(key = paste(Line, CHROM, POS1, sep = "_")) %>%
    select(Line, CHROM, CHROM_label, POS1, key)
  
  shared_keys <- intersect(early_df$key, late_df$key)
  
  early_only <- early_df %>% filter(!key %in% shared_keys) %>% mutate(Type = "Early only")
  late_only  <- late_df  %>% filter(!key %in% shared_keys) %>% mutate(Type = "Late only")
  shared_pts <- early_df %>% filter( key %in% shared_keys) %>% mutate(Type = "Shared")
  
  plot_df <- bind_rows(early_only, late_only, shared_pts)
  plot_df$Type <- factor(plot_df$Type,
                         levels = c("Early only", "Late only", "Shared"))
  
  # Plot crossovers
  p_cross <- ggplot(plot_df, aes(x = POS1 / 1e6, y = Line, color = Type)) +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(~ CHROM_label, ncol = 4, scales = "free_x") +
    scale_color_manual(
      values = c("Early only" = "#FF7F7F",
                 "Late only"  = "skyblue",
                 "Shared"     = "purple"),
      name = "Crossover breakpoint"
    ) +
    labs(x = "Position within chromosome (Mb)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
  
  # Summary tables
  overall_summary <- data.frame(
    Category = c("Early only", "Late only", "Shared (early + late)"),
    n_events = c(nrow(early_only), nrow(late_only), nrow(shared_pts))
  )
  overall_summary$percent_of_total <- round(
    overall_summary$n_events / sum(overall_summary$n_events) * 100, 2)
  
  overlap_rate_vs_early <- nrow(shared_pts) /
    (nrow(early_only) + nrow(shared_pts)) * 100
  
  overall_summary <- rbind(
    overall_summary,
    data.frame(
      Category = "Overlap rate vs early breakpoints (%)",
      n_events = round(overlap_rate_vs_early, 2),
      percent_of_total = NA
    )
  )
  
  chrom_summary <- plot_df %>%
    group_by(CHROM_label, Type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Type, values_from = n, values_fill = 0) %>%
    mutate(
      Total = `Early only` + `Late only` + Shared,
      Overlap_percent = round(Shared / Total * 100, 2)
    ) %>%
    arrange(CHROM_label)
  
  line_summary <- plot_df %>%
    group_by(Line, Type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Type, values_from = n, values_fill = 0) %>%
    mutate(
      Total = `Early only` + `Late only` + Shared,
      Overlap_percent = round(Shared / Total * 100, 2)
    ) %>%
    arrange(Line)
  
  print(overall_summary)
  print(chrom_summary)
  print(line_summary)
  
  write.table(overall_summary, file.path(output_dir, "SOLine_timecourse_crossover_overlap_overall.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(chrom_summary, file.path(output_dir, "SOLine_timecourse_crossover_overlap_by_chrom.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(line_summary, file.path(output_dir, "SOLine_timecourse_crossover_overlap_by_line.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Save plots and summary
ggsave(file.path(output_dir, "SOLine_timecourse_delta_hist.pdf"), p_hist, width = 14, height = 8)
ggsave(file.path(output_dir, "SOLine_timecourse_correlation.pdf"), p_corr, width = 10, height = 5)
ggsave(file.path(output_dir, "SOLine_timecourse_stable_fraction.pdf"), p_stable, width = 10, height = 5)
if (exists("p_cross")) {
  ggsave(file.path(output_dir, "SOLine_timecourse_crossover_compare.pdf"), p_cross, width = 14, height = 8)
}

write.table(summary_stats,
            file = file.path(output_dir, "SOLine_timecourse_stability_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("SOLine time-course haplotype stability analysis completed.\n")