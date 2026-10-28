# ============================================================
# Heterogeneity and Coverage Plot with Three Depth Levels
# ============================================================

library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

# ============================================================
# 1. Load data
# ============================================================

df <- read.table(file = "/data/01_Het_Fws/het-cov.txt", header = T, row.names = NULL)

# ============================================================
# 2. Data preprocessing
# ============================================================

# Reshape data for coverage levels
df_long <- df %>%
  select(ID, Heterozygosity, Coverage_1x, Coverage_10x, Coverage_30x, Category) %>%
  pivot_longer(
    cols = c(Coverage_1x, Coverage_10x, Coverage_30x),
    names_to = "Coverage_Level",
    values_to = "Coverage_Value"
  ) %>%
  mutate(
    Coverage_Level = factor(Coverage_Level, 
                            levels = c("Coverage_1x", "Coverage_10x", "Coverage_30x"),
                            labels = c("1x", "10x", "30x"))
  )

# ============================================================
# 3. Color palette
# ============================================================

category_colors <- c("Single_oocyst" = "#FFD456", "Single_sporozoite" = "#7F88CB")

coverage_colors <- c("1x" = "#E57373", "10x" = "#E57373", "30x" = "#E57373")

coverage_linetypes <- c("1x" = "dotted", "10x" = "dashed", "30x" = "solid")

# ============================================================
# 4. Create plot
# ============================================================

p <- ggplot() +
  # 1. Bar plot for Heterozygosity
  geom_col(data = df, 
           aes(x = ID, y = Heterozygosity, fill = Category), 
           width = 0.7, alpha = 0.85) +
  
  # 2. Lines and points for three coverage levels (color and linetype both mapped)
  geom_line(data = df_long, 
            aes(x = ID, y = Coverage_Value / 100, 
                color = Coverage_Level, linetype = Coverage_Level, group = Coverage_Level),
            linewidth = 0.8) +
  geom_point(data = df_long, 
             aes(x = ID, y = Coverage_Value / 100, 
                 color = Coverage_Level, shape = Coverage_Level),  
             size = 2.5) +
  
  # 3. Dual y-axis
  scale_y_continuous(
    name = "Heterogeneity",
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.2),
    sec.axis = sec_axis(
      ~ . * 100,
      name = "Breadth of coverage (%)",
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%")
    )
  ) +
  
  # 4. Color and linetype settings
  scale_fill_manual(values = category_colors, name = "Category") +
  scale_color_manual(values = coverage_colors, name = "Coverage Depth") +   
  scale_linetype_manual(values = coverage_linetypes, name = "Coverage Depth") +  
  scale_shape_manual(values = c("1x" = 16, "10x" = 17, "30x" = 18), name = "Coverage Depth") +
  
  # 5. Theme
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y.left = element_text(color = "black", size = 11),
    axis.title.y.left = element_text(color = "black", face = "bold", size = 12),
    axis.text.y.right = element_text(color = "#E57373", size = 11),     
    axis.title.y.right = element_text(color = "#E57373", face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = unit(c(5, 5, 5, 5), "mm")
  ) +
  labs(
    title = "",
    x = "Sample ID",
    fill = "Category",
    color = "Coverage Depth",
    linetype = "Coverage Depth",
    shape = "Coverage Depth"
  )

# ============================================================
# 5. Display and save
# ============================================================

print(p)

ggsave("/results/01_Het_Fws/het_coverage_plot/heterozygosity_coverage_plot.png", plot = p, width = 10, height = 6, dpi = 300)
ggsave("/results/01_Het_Fws/het_coverage_plot/heterozygosity_coverage_plot.pdf", plot = p, width = 10, height = 6)
