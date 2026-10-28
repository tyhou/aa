# Recombination map (physical vs genetic distance)
# Inputs: /data/05_Genetic_Map
# Outputs: /results/05_Genetic_Map

library(dplyr)
library(ggplot2)

# Paths
input_dir  <- "/data/05_Genetic_Map"
output_dir <- "/results/05_Genetic_Map"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
rho <- read.table(file.path(input_dir, "Rho.stat.final.txt"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(rho) <- c("Chr", "Loci", "Mean_rho", "Median", "L95", "U95")

# Clean
rho <- rho %>% filter(Loci > 0)
rho$Loci <- as.numeric(rho$Loci)
rho$Median <- as.numeric(rho$Median)
rho$Chr <- factor(rho$Chr, levels = paste0("Chr", 1:8))
rho <- rho %>% arrange(Chr, Loci)

# Cumulative distances
cumulative_median <- rho %>%
  group_by(Chr) %>%
  mutate(
    physical_cum = Loci - first(Loci),
    delta_phys   = c(0, diff(Loci)),
    genetic_cum  = cumsum(Median * delta_phys)
  ) %>%
  ungroup() %>%
  select(Chr, Loci, Median, physical_cum, genetic_cum)

write.table(cumulative_median, file.path(output_dir, "cumulative_physical_genetic_median.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plotting prep
cum <- cumulative_median %>%
  mutate(physical_cum = as.numeric(physical_cum), genetic_cum = as.numeric(genetic_cum)) %>%
  arrange(Chr, physical_cum) %>%
  filter(!is.na(physical_cum), !is.na(genetic_cum))

max_phys <- max(cum$physical_cum, na.rm = TRUE)
max_gen  <- max(cum$genetic_cum,  na.rm = TRUE)
scale_factor <- max_phys / max_gen

cum <- cum %>% mutate(genetic_scaled = genetic_cum * scale_factor)

chr_max <- cum %>%
  group_by(Chr) %>%
  summarise(max_phys = max(physical_cum, na.rm = TRUE), max_gen = max(genetic_cum, na.rm = TRUE)) %>%
  ungroup()

chr_bar <- bind_rows(
  chr_max %>% mutate(group = "physical", y = max_phys),
  chr_max %>% mutate(group = "genetic",  y = max_gen * scale_factor)
)

seg_df <- cum %>%
  group_by(Chr) %>%
  mutate(gap = c(0, diff(physical_cum))) %>%
  ungroup() %>%
  filter(gap > 0.5) %>%
  select(Chr, physical_cum, genetic_scaled) %>%
  mutate(x = "physical", xend = "genetic", y = physical_cum, yend = genetic_scaled)

phys_breaks <- seq(0, ceiling(max_phys / 500) * 500, by = 500)
gen_breaks_raw <- seq(0, ceiling(max_gen / 100) * 100, by = 100)
gen_breaks_scaled <- gen_breaks_raw * scale_factor
y_max <- max(max_phys, max(gen_breaks_scaled)) * 1.05

g <- ggplot() +
  geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend, group = Chr),
               color = "gray70", linewidth = 0.35, alpha = 0.7) +
  geom_col(data = chr_bar, aes(x = group, y = y, fill = group), position = "dodge", width = 0.12) +
  geom_line(data = cum, aes(x = "physical", y = physical_cum, group = Chr), color = "#3949AB", linewidth = 2) +
  geom_line(data = cum, aes(x = "genetic", y = genetic_scaled, group = Chr), color = "#B71C1C", linewidth = 2) +
  facet_wrap(~ Chr, strip.position = "bottom", ncol = 8) +
  labs(x = "", y = "Genetic distance (cM)") +
  scale_y_continuous(
    limits = c(0, y_max), breaks = gen_breaks_scaled, labels = gen_breaks_raw,
    expand = c(0.02, 0.02),
    sec.axis = sec_axis(~ ., name = "Physical distance (kb)", breaks = phys_breaks)
  ) +
  scale_fill_manual(values = c("physical" = "#3949AB", "genetic" = "#B71C1C")) +
  theme(
    panel.border = element_blank(), strip.background = element_blank(),
    axis.ticks.x = element_blank(), panel.grid = element_blank(),
    panel.background = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12), strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "cumulative_physical_genetic.pdf"), g, width = 7, height = 4, dpi = 300)
ggsave(file.path(output_dir, "cumulative_physical_genetic.png"), g, width = 7, height = 4, dpi = 600)

cat("Recombination map analysis completed.\n")