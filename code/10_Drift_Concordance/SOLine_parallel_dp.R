# Paths adapted for this repository (see README): inputs from /data, outputs to /results.
input_dir  <- "/data/04_Allele_Frequency/SOLine"
gff_dir    <- "/data/06_Fst_Pi"
output_dir <- "/results/10_Drift_Concordance"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
# 安装需要的包（如果已经安装可以跳过）
# install.packages(c("ggplot2", "dplyr", "tidyr"))

library(ggplot2)
library(dplyr)
library(tidyr)

# ========== 1. 读取数据 ==========
freq <- read.table(file.path(input_dir, "SOLine_alt_allele_frequency_1391.tsv"),
                   header = TRUE,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)

# ========== 2. 定义早晚样品对 ==========
line_pairs <- data.frame(
  Line  = c("Line1","Line2","Line3","Line4","Line6","Line7","Line8","Line10","Line11","Line16","Line19","Line24","Line30","Line31","Line32"),
  Early = c("W683","W685","W373","W688","W469","W376","W447","W377","W378","W379","W380","W441","W388","W383","W384"),
  Late  = c("W455","W454","W468","W689","W375","W456","W453","W457","W470","W458","W459","W446","W464","W465","W466")
)

# 检查这些列是否都存在
all_samples <- c(line_pairs$Early, line_pairs$Late)
missing <- setdiff(all_samples, colnames(freq))
if (length(missing) > 0) {
  stop(paste("以下样品列在文件中不存在：", paste(missing, collapse = ", ")))
}

# ========== 3. 提取需要的列 ==========
freq_sub <- freq[, c("CHROM", "POS", "ID", "REF", "ALT", all_samples)]

# ========== 4. 计算每个 SNP 在每个 Line 的 Δp ==========
results <- data.frame()

for (i in 1:nrow(line_pairs)) {
  line  <- line_pairs$Line[i]
  early <- line_pairs$Early[i]
  late  <- line_pairs$Late[i]
  
  dp <- freq_sub[[late]] - freq_sub[[early]]
  
  temp <- data.frame(
    CHROM = freq_sub$CHROM,
    POS   = freq_sub$POS,
    ID    = freq_sub$ID,
    REF   = freq_sub$REF,
    ALT   = freq_sub$ALT,
    Line  = line,
    dp    = dp
  )
  results <- rbind(results, temp)
}

# ========== 5. 汇总每个 SNP 的平行性 ==========
summary_res <- results %>%
  group_by(CHROM, POS, ID, REF, ALT) %>%
  summarise(
    n         = n(),
    n_up      = sum(dp > 0, na.rm = TRUE),
    n_down    = sum(dp < 0, na.rm = TRUE),
    mean_dp   = mean(dp, na.rm = TRUE),
    median_dp = median(dp, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    n_same      = pmax(n_up, n_down),
    same_frac   = n_same / n,
    abs_mean_dp = abs(mean_dp)
  ) %>%
  arrange(desc(same_frac), desc(abs_mean_dp))

# 输出汇总表格
write.csv(summary_res, "parallel_dp_summary.csv", row.names = FALSE)

# 按 0.2 分档
summary_res <- summary_res %>%
  mutate(bin = cut(same_frac,
                   breaks = seq(0, 1, by = 0.2),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")))

# 统计每档位点数
bin_count <- summary_res %>%
  group_by(bin) %>%
  summarise(n = n()) %>%
  mutate(pct = n / sum(n) * 100)

print(bin_count)

# 画图：纵坐标改为比例
p <- ggplot(bin_count, aes(x = bin, y = pct)) +
  geom_col(fill = "grey60", color = "black") +
  geom_text(aes(label = paste0(round(pct, 1))),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(title = "Distribution of cross-line direction concordance",
       x = "same_frac", y = "Proportion of SNPs (%)") +
  theme_minimal()

ggsave("same_frac_bin_pct.pdf", p, width = 6, height = 5)
print(p)

# ========== 6. 候选位点基因注释（不用 rtracklayer） ==========

# 筛选 same_frac >= 0.8 的候选位点
candidate_snps <- summary_res %>% filter(same_frac >= 0.8)
cat("候选位点数：", nrow(candidate_snps), "\n")

# ---- 读取 GFF ----
gff_file <- file.path(gff_dir, "IId_11730_genome.gff3")

gff <- read.delim(gff_file, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, quote = "",
                  col.names = c("seqid","source","type","start","end",
                                "score","strand","phase","attributes"))

cat("GFF 中的 type：\n"); print(unique(gff$type))
cat("GFF 中的染色体：\n"); print(unique(gff$seqid))

# ---- 选择基因注释 ----
genes <- gff[gff$type == "gene", ]
if (nrow(genes) == 0) genes <- gff[gff$type == "mRNA", ]
cat("用于注释的基因数：", nrow(genes), "\n")

# ---- 从 attributes 列提取基因 ID / locus_tag / biotype ----
# 修复版：始终返回与输入等长的向量，找不到匹配则返回 NA
extract_attr <- function(attr_vec, key) {
  pattern <- paste0("(^|;)\\s*", key, "\\s*=\\s*([^;]*)")
  m <- regexec(pattern, attr_vec)
  res <- regmatches(attr_vec, m)
  sapply(res, function(x) {
    if (length(x) >= 3) return(x[3]) else return(NA_character_)
  })
}

genes$gene_id   <- extract_attr(genes$attributes, "ID")
genes$gene_name <- extract_attr(genes$attributes, "locus_tag")
genes$biotype   <- extract_attr(genes$attributes, "gene_biotype")

# 若 locus_tag 缺失，用 ID 去掉 "gene-" 前缀作为备份
genes$gene_name[is.na(genes$gene_name)] <-
  sub("^gene-", "", genes$gene_id[is.na(genes$gene_name)])

# 检查一下
cat("基因属性示例：\n")
print(head(genes[, c("seqid","start","end","strand","gene_id","gene_name","biotype")], 5))

# ---- 染色体命名匹配 ----
match_chrom <- function(chrom, gff_chroms) {
  if (chrom %in% gff_chroms) return(chrom)
  if (paste0("Chr", chrom) %in% gff_chroms) return(paste0("Chr", chrom))
  if (paste0("chr", chrom) %in% gff_chroms) return(paste0("chr", chrom))
  if (paste0("chr", chrom, "_") %in% gff_chroms) return(paste0("chr", chrom, "_"))
  if (grepl("^Chr", chrom) && sub("^Chr", "", chrom) %in% gff_chroms) return(sub("^Chr", "", chrom))
  if (grepl("^chr", chrom) && sub("^chr", "", chrom) %in% gff_chroms) return(sub("^chr", "", chrom))
  return(NA)
}

gff_chroms <- unique(gff$seqid)
candidate_snps$CHROM_gff <- sapply(as.character(candidate_snps$CHROM),
                                   match_chrom, gff_chroms = gff_chroms)
cat("匹配失败的染色体：\n")
print(unique(candidate_snps$CHROM[is.na(candidate_snps$CHROM_gff)]))

# ---- 逐个 SNP 查找包含它的基因 ----
annot_list <- list()

for (i in seq_len(nrow(candidate_snps))) {
  chrom <- candidate_snps$CHROM_gff[i]
  pos   <- candidate_snps$POS[i]
  
  if (is.na(chrom)) {
    annot_list[[i]] <- data.frame(
      CHROM = candidate_snps$CHROM[i], POS = pos,
      gene_id = NA, gene_name = NA, biotype = NA,
      gene_start = NA, gene_end = NA, strand = NA,
      distance_to_gene = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  sub_genes <- genes[genes$seqid == chrom, ]
  
  # 先找包含该位点的基因
  hit <- sub_genes[sub_genes$start <= pos & sub_genes$end >= pos, ]
  
  if (nrow(hit) > 0) {
    annot_list[[i]] <- data.frame(
      CHROM = candidate_snps$CHROM[i], POS = pos,
      gene_id    = hit$gene_id[1],
      gene_name  = hit$gene_name[1],
      biotype    = hit$biotype[1],
      gene_start = hit$start[1],
      gene_end   = hit$end[1],
      strand     = hit$strand[1],
      distance_to_gene = 0,
      stringsAsFactors = FALSE
    )
  } else {
    # 位点不在任何基因内，计算到最近基因的距离
    if (nrow(sub_genes) > 0) {
      dist_up   <- pos - sub_genes$end
      dist_down <- sub_genes$start - pos
      dist_all  <- ifelse(dist_up >= 0, dist_up,
                          ifelse(dist_down >= 0, dist_down, 0))
      nearest   <- which.min(abs(dist_all))
      annot_list[[i]] <- data.frame(
        CHROM = candidate_snps$CHROM[i], POS = pos,
        gene_id    = sub_genes$gene_id[nearest],
        gene_name  = sub_genes$gene_name[nearest],
        biotype    = sub_genes$biotype[nearest],
        gene_start = sub_genes$start[nearest],
        gene_end   = sub_genes$end[nearest],
        strand     = sub_genes$strand[nearest],
        distance_to_gene = dist_all[nearest],
        stringsAsFactors = FALSE
      )
    } else {
      annot_list[[i]] <- data.frame(
        CHROM = candidate_snps$CHROM[i], POS = pos,
        gene_id = NA, gene_name = NA, biotype = NA,
        gene_start = NA, gene_end = NA, strand = NA,
        distance_to_gene = NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

annot <- do.call(rbind, annot_list)

# ---- 合并回候选位点表 ----
annot <- merge(candidate_snps, annot, by = c("CHROM","POS"), all.x = TRUE)

# ---- 输出 ----
write.csv(annot, "candidate_snp_genes.csv", row.names = FALSE)
cat("注释结果已保存至 candidate_snp_genes.csv\n")
print(annot)
