# Paths adapted for this repository (see README): inputs from /data, outputs to /results.
input_dir  <- "/data/11_AlphaFold"
output_dir <- "/results/11_AlphaFold"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
# ============================================
# 1. 加载必要的包（都已安装过）
# ============================================
library(jsonlite)
library(ggplot2)
library(data.table)
library(scales)
library(RColorBrewer)

# ============================================
# 2. 读取并提取数据
# ============================================
data <- fromJSON(file.path(input_dir, "fold_2024_12_02_21_37_cgd5_4090_hlj_full_data_0.json"),
                 simplifyVector = TRUE)

print(names(data))  # 确认有 pae, token_res_ids, token_chain_ids

pae_mat     <- data$pae
token_res   <- data$token_res_ids
token_chain <- data$token_chain_ids

cat("PAE 矩阵维度:", dim(pae_mat), "\n")
cat("Token 数量:", length(token_res), "\n")

# 保存原始数据（作为补充材料）
write.csv(pae_mat, "pae_matrix.csv", row.names = FALSE)
token_map <- data.frame(
  token_index = 1:length(token_res),
  chain       = token_chain,
  resid       = token_res
)
write.csv(token_map, "token_mapping.csv", row.names = FALSE)

# ============================================
# 3. 定义结构域边界（来自 SMART/PFAM 截图）
# ============================================
domain_dexdc  <- 261:479
domain_helicc <- 682:777
domain_ha2    <- 838:942

# 找到各结构域对应的 token 索引
idx_dexdc  <- which(token_res %in% domain_dexdc)
idx_helicc <- which(token_res %in% domain_helicc)
idx_ha2    <- which(token_res %in% domain_ha2)

# ============================================
# 4. 计算域内 / 域间平均 PAE
# ============================================
pae_dexdc_intra  <- mean(pae_mat[idx_dexdc,  idx_dexdc],  na.rm = TRUE)
pae_helicc_intra <- mean(pae_mat[idx_helicc, idx_helicc], na.rm = TRUE)
pae_ha2_intra    <- mean(pae_mat[idx_ha2,    idx_ha2],    na.rm = TRUE)

pae_dexdc_helicc <- mean(pae_mat[idx_dexdc,  idx_helicc], na.rm = TRUE)
pae_dexdc_ha2    <- mean(pae_mat[idx_dexdc,  idx_ha2],    na.rm = TRUE)
pae_helicc_ha2   <- mean(pae_mat[idx_helicc, idx_ha2],    na.rm = TRUE)

cat("\n=== PAE 统计结果 ===\n")
cat("DEXDc (261-479) 域内 PAE:",  round(pae_dexdc_intra,  2), "Å\n")
cat("HELICc (682-777) 域内 PAE:", round(pae_helicc_intra, 2), "Å\n")
cat("HA2 (838-942) 域内 PAE:",    round(pae_ha2_intra,    2), "Å\n")
cat("DEXDc-HELICc 域间 PAE:",     round(pae_dexdc_helicc, 2), "Å\n")
cat("DEXDc-HA2 域间 PAE:",        round(pae_dexdc_ha2,    2), "Å\n")
cat("HELICc-HA2 域间 PAE:",       round(pae_helicc_ha2,   2), "Å\n")

pae_stats <- data.frame(
  Region   = c("DEXDc intra", "HELICc intra", "HA2 intra",
               "DEXDc-HELICc inter", "DEXDc-HA2 inter", "HELICc-HA2 inter"),
  Mean_PAE = c(pae_dexdc_intra, pae_helicc_intra, pae_ha2_intra,
               pae_dexdc_helicc, pae_dexdc_ha2, pae_helicc_ha2)
)
write.csv(pae_stats, "pae_domain_stats.csv", row.names = FALSE)

# ============================================
# 5. 用 ggplot2 画 PAE 热图（带结构域方框）
# ============================================

# 把 PAE 矩阵转为长格式
# pae_mat[i, j] 对应 row = i, col = j
# as.vector(pae_mat) 是按列展开：j=1 时 i=1..n，j=2 时 i=1..n ...
n <- nrow(pae_mat)
pae_long <- data.table(
  x   = rep(1:n, times = n),  # col 索引（Aligned residue）
  y   = rep(1:n, each  = n),  # row 索引（Scored residue）
  pae = as.vector(pae_mat)
)

# 结构域方框
domain_rects <- data.frame(
  xmin  = c(min(idx_dexdc), min(idx_helicc), min(idx_ha2)) - 0.5,
  xmax  = c(max(idx_dexdc), max(idx_helicc), max(idx_ha2)) + 0.5,
  ymin  = c(min(idx_dexdc), min(idx_helicc), min(idx_ha2)) - 0.5,
  ymax  = c(max(idx_dexdc), max(idx_helicc), max(idx_ha2)) + 0.5,
  label = c("DEXDc", "HELICc", "HA2")
)

p <- ggplot(pae_long, aes(x = x, y = y, fill = pae)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = c("darkgreen", "#80d35d", "#ddfad1"),
    values = scales::rescale(c(0, 10, 20, 30)),
    limits = c(0, 30),
    oob    = scales::squish,
    name   = "PAE (Å)"
  ) +
  # 画结构域方框
  geom_rect(
    data = domain_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, color = "white", fill = NA, linewidth = 0.5
  ) +
  # 在方框中心加文字标签
  geom_text(
    data = domain_rects,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    inherit.aes = FALSE, color = "white", size = 1.8, fontface = "bold"
  ) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_fixed() +
  labs(
    x = "Aligned residue",
    y = "Scored residue",
    title = "PAE matrix of CpDHR1"
  ) +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    plot.title  = element_text(hjust = 0.5),
    axis.text   = element_text(size = 9)
  )

# 保存为 PNG
ggsave("pae_heatmap_with_domains.png", p,
       width = 5, height = 5, dpi = 300)

# 保存为 PDF
ggsave("pae_heatmap_with_domains.pdf", p,
       width = 5, height = 5)

cat("\n绘图完成：pae_heatmap_with_domains.png\n")


