# Paths adapted for this repository (see README): inputs from /data, outputs to /results.
input_dir  <- "/data/11_AlphaFold"
output_dir <- "/results/11_AlphaFold"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
library(ggplot2)
library(scales)

# 创建一个 0 到 100 的 pLDDT 数据
df_legend <- data.frame(
  x = 1,
  y = seq(0, 100, length.out = 100),
  plddt = seq(0, 100, length.out = 100)
)

# 绘制渐变图例
p_legend <- ggplot(df_legend, aes(x = x, y = y, fill = plddt)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#FF7D45", "#FFDB13", "#65CBF3", "#0053D6"), # 官方四色
    values = scales::rescale(c(0, 50, 70, 90, 100)),        # 刻度映射
    limits = c(0, 100),
    name = "pLDDT"
  ) +
  scale_y_continuous(
    breaks = c(0, 50, 70, 90, 100), # 标出关键阈值
    position = "right"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() + # 去掉所有背景和坐标轴
  theme(
    legend.position = "none", # 我们画的就是图例本身，不需要额外的图例
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# 保存为透明背景的 PNG，方便后期排版
ggsave("pLDDT_legend.png", p_legend, width = 1.5, height = 4, dpi = 300, bg = "transparent")
ggsave("pLDDT_legend.pdf", p_legend, width = 1.5, height = 4, bg = "transparent")

