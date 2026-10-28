# Paths adapted for this repository (see README): inputs from /data, outputs to /results.
input_dir  <- "/data/11_AlphaFold"
output_dir <- "/results/11_AlphaFold"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
# ==============================================================================
# AlphaFold 3 pLDDT 分析脚本（不依赖 bio3d）
# ==============================================================================

# 1. 加载R包
library(jsonlite)
library(ggplot2)
library(dplyr)
library(scales)

# ==============================================================================
# 2. 文件路径和结构域定义
# ==============================================================================

json_file <- file.path(input_dir, "fold_2024_12_02_21_37_cgd5_4090_hlj_full_data_0.json")
cif_file  <- file.path(input_dir, "fold_2024_12_02_21_37_cgd5_4090_hlj_model_0.cif")

# 定义结构域边界（根据SMART/PFAM分析）
domain_definitions <- list(
  "N-terminal" = c(1, 260),
  "DEXDc"      = c(261, 479),
  "Linker1"    = c(480, 681),
  "HELICc"     = c(682, 777),
  "Linker2"    = c(778, 837),
  "HA2"        = c(838, 942),
  "C-terminal" = c(943, 1272)  # 请根据实际序列长度调整
)

# ==============================================================================
# 3. 从CIF文件解析原子记录（不依赖bio3d）
# ==============================================================================

cat("正在解析CIF文件...\n")

# 读取所有行
cif_lines <- readLines(cif_file, warn = FALSE)

# 找到 _atom_site 循环的起始位置
atom_site_start <- grep("^_atom_site\\.", cif_lines)
if (length(atom_site_start) == 0) stop("未在CIF文件中找到 _atom_site 记录")

# 提取列名（去掉 _atom_site. 前缀）
col_names <- sub("^_atom_site\\.", "", cif_lines[atom_site_start])
# 找数据行的起始位置（列名之后的第一个非注释行）
data_start <- max(atom_site_start) + 1
while (data_start <= length(cif_lines) && grepl("^#|^_|^loop_|^$", cif_lines[data_start])) {
  data_start <- data_start + 1
}
# 找数据结束位置（遇到 # 或空行或下一个 loop_）
data_end <- data_start
while (data_end <= length(cif_lines) && !grepl("^#|^loop_|^_", cif_lines[data_end]) && nzchar(trimws(cif_lines[data_end]))) {
  data_end <- data_end + 1
}
data_end <- data_end - 1

# 提取原子数据
atom_lines <- cif_lines[data_start:data_end]
# 按空白分割每一行
atom_list <- strsplit(trimws(atom_lines), "\\s+")

# 找出需要的列索引
col_idx_resid <- which(col_names == "label_seq_id")
col_idx_chain <- which(col_names == "label_asym_id")
col_idx_bfac  <- which(col_names == "B_iso_or_equiv")
col_idx_atom  <- which(col_names == "label_atom_id")

if (length(col_idx_bfac) == 0) stop("未在CIF文件中找到 B_iso_or_equiv 列")

# 构建数据框
atom_info <- data.frame(
  chain = sapply(atom_list, `[`, col_idx_chain),
  resid = as.integer(sapply(atom_list, `[`, col_idx_resid)),
  bfactor = as.numeric(sapply(atom_list, `[`, col_idx_bfac)),
  stringsAsFactors = FALSE
)

cat("解析完成，共", nrow(atom_info), "个原子\n")

# ==============================================================================
# 4. 按残基计算平均pLDDT
# ==============================================================================

cat("正在计算逐残基pLDDT...\n")
res_plddt <- atom_info %>%
  group_by(chain, resid) %>%
  summarise(plddt_mean = mean(bfactor, na.rm = TRUE), .groups = "drop") %>%
  arrange(chain, resid)

# 保存逐残基pLDDT为CSV（补充数据）
write.csv(res_plddt, "Supplementary_Data_pLDDT_per_residue.csv", row.names = FALSE)

# ==============================================================================
# 5. 绘制逐残基pLDDT曲线图
# ==============================================================================

cat("正在绘制逐残基pLDDT曲线...\n")

# 添加结构域注释
res_plddt <- res_plddt %>%
  mutate(domain = case_when(
    resid >= domain_definitions$`N-terminal`[1] & resid <= domain_definitions$`N-terminal`[2] ~ "N-terminal",
    resid >= domain_definitions$DEXDc[1] & resid <= domain_definitions$DEXDc[2] ~ "DEXDc",
    resid >= domain_definitions$Linker1[1] & resid <= domain_definitions$Linker1[2] ~ "Linker1",
    resid >= domain_definitions$HELICc[1] & resid <= domain_definitions$HELICc[2] ~ "HELICc",
    resid >= domain_definitions$Linker2[1] & resid <= domain_definitions$Linker2[2] ~ "Linker2",
    resid >= domain_definitions$HA2[1] & resid <= domain_definitions$HA2[2] ~ "HA2",
    resid >= domain_definitions$`C-terminal`[1] & resid <= domain_definitions$`C-terminal`[2] ~ "C-terminal",
    TRUE ~ "Other"
  ))

domain_colors <- c(
  "N-terminal" = "grey50", "DEXDc" = "#008787", "Linker1" = "grey50",
  "HELICc" = "#ACD7D7", "Linker2" = "grey50", "HA2" = "#2D882D",
  "C-terminal" = "grey50"
)

# 计算每个结构域的xmin和xmax用于画背景色块
domain_rects <- res_plddt %>%
  group_by(domain) %>%
  summarise(xmin = min(resid), xmax = max(resid), .groups = "drop")

p_plddt <- ggplot(res_plddt, aes(x = resid, y = plddt_mean)) +
  geom_rect(data = domain_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = domain),
            inherit.aes = FALSE, alpha = 0.15) +
  geom_line(color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(50, 70, 90), linetype = "dashed",
             color = c("orange", "yellow", "blue"), alpha = 0.7) +
  scale_fill_manual(values = domain_colors, name = "Domain") +
  labs(x = "Residue number", y = "pLDDT",
       title = "Per-residue pLDDT of CpDHR1 (AlphaFold 3)") +
  ylim(0, 100) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("Supplementary_Figure_pLDDT_profile.png", p_plddt, width = 12, height = 5, dpi = 300)
ggsave("Supplementary_Figure_pLDDT_profile.pdf", p_plddt, width = 12, height = 5)

# ==============================================================================
# 6. 计算各结构域的平均pLDDT并制表
# ==============================================================================

cat("正在计算各结构域的平均pLDDT...\n")

domain_stats <- res_plddt %>%
  group_by(domain) %>%
  summarise(
    Residue_Range = paste0(min(resid), "-", max(resid)),
    Mean_pLDDT    = round(mean(plddt_mean, na.rm = TRUE), 1),
    Median_pLDDT  = round(median(plddt_mean, na.rm = TRUE), 1),
    Min_pLDDT     = round(min(plddt_mean, na.rm = TRUE), 1),
    Max_pLDDT     = round(max(plddt_mean, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(match(domain, names(domain_definitions)))

write.csv(domain_stats, "Supplementary_Table_pLDDT_by_domain.csv", row.names = FALSE)

cat("\n分析完成！已生成：\n")
cat("1. Supplementary_Data_pLDDT_per_residue.csv\n")
cat("2. Supplementary_Figure_pLDDT_profile.png/pdf\n")
cat("3. Supplementary_Table_pLDDT_by_domain.csv\n")

print(domain_stats)

