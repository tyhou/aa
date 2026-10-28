# Paths adapted for this repository (see README): inputs from /data, outputs to /results.
input_dir  <- "/data/04_Allele_Frequency/SOLine"
output_dir <- "/results/01_Het_Fws/Het_dotplot"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
library(readxl)
library(ggplot2)


# 读取数据
df <- read_excel(file.path(input_dir, "Het-SOLine-input.xlsx"))

# 关键：去掉 Line 列中的空格，使 "Line 1" 变为 "Line1" 以匹配颜色映射
df$Line <- gsub(" ", "", df$Line)

# 颜色映射（已去掉 HLJ, GD, F1）
group_order <- c("Line1", "Line2", "Line3", "Line4", "Line6",
                 "Line7", "Line8", "Line10", "Line11", "Line16",
                 "Line19", "Line24", "Line30", "Line31", "Line32")
group_colors <- setNames(c(
  "#9BC2E6", "#F4B183", "#C15811",
  "#757575", "#F06292", "#BA124A", "#BA8CDC", "#7030A0", "#BF7F3F",
  "#0964AF", "#774F27", "#063F6E", "#68A141", "#02192C", "#385723"
), group_order)

# Line 设为因子，固定顺序
df$Line <- factor(df$Line, levels = group_order)

# 横坐标：每 3 一个刻度线，但只在 6 的倍数显示标签
xbreaks <- seq(0, 60, 3)
xlabels <- ifelse(xbreaks %% 6 == 0, as.character(xbreaks), "")

# 散点图
p <- ggplot(df, aes(x = DPI, y = Het, color = Line)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = group_colors, breaks = group_order) +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Days post infection", y = "Heterozygosity") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text  = element_text(color = "black", size = 13),  # 刻度标签黑色，与标题同大小
    legend.key  = element_rect(fill = "white", colour = NA),
    legend.text = element_text(size = 10)
  )
p

# 保存
ggsave("Het-SOLine-dotplot.png", p, width = 6, height = 3, dpi = 300)
ggsave("Het-SOLine-dotplot.pdf", p, width = 6, height = 3)

print(p)
