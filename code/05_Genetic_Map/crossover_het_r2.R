# Crossover vs Heterozygosity correlation
# Inputs: /data/05_Genetic_Map
# Outputs: /results/05_Genetic_Map

library(dplyr)
library(ggplot2)

# Paths
input_dir  <- "/data/05_Genetic_Map"
output_dir <- "/results/05_Genetic_Map"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
merged_ok <- read.csv(file.path(input_dir, "crossover_vs_het_data.csv"), header = TRUE, stringsAsFactors = FALSE)

if (!"Total_Crossover" %in% colnames(merged_ok) || !"Het_Rate" %in% colnames(merged_ok)) {
  stop("crossover_vs_het_data.csv must contain 'Total_Crossover' and 'Het_Rate' columns")
}

merged_ok <- merged_ok %>% filter(!is.na(Total_Crossover), !is.na(Het_Rate))

# Correlation
pearson  <- cor.test(merged_ok$Total_Crossover, merged_ok$Het_Rate, method = "pearson")
r2_pearson <- pearson$estimate^2

cat("Pearson r =", round(pearson$estimate, 4), "\n")
cat("R2 =", round(r2_pearson, 4), "\n")
cat("p =", format.pval(pearson$p.value, digits = 4), "\n")

# Plot
p_label <- paste0("p = ", signif(pearson$p.value, 3))
r2_label <- paste0("R\u00B2 = ", round(r2_pearson, 3))
label_text <- paste0(r2_label, "\n", p_label)

max_x <- max(merged_ok$Total_Crossover, na.rm = TRUE)
max_y <- max(merged_ok$Het_Rate, na.rm = TRUE)

p <- ggplot(merged_ok, aes(x = Total_Crossover, y = Het_Rate)) +
  geom_point(size = 3, color = "#3949AB", alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#3949AB", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max_x * 0.95, y = max_y * 0.95, label = label_text, hjust = 1, vjust = 1, size = 4.5) +
  coord_cartesian(xlim = c(0, max_x * 1.05), ylim = c(0, max_y * 1.05)) +
  scale_x_continuous(breaks = pretty(merged_ok$Total_Crossover, n = 5)) +
  scale_y_continuous(breaks = pretty(merged_ok$Het_Rate, n = 5)) +
  labs(x = "Total Crossover Count", y = "Heterozygosity") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

ggsave(file.path(output_dir, "crossover_vs_het.png"), p, width = 5, height = 4, dpi = 300)
ggsave(file.path(output_dir, "crossover_vs_het.pdf"), p, width = 5, height = 4, device = "pdf")

cat("Crossover vs Heterozygosity analysis completed.\n")