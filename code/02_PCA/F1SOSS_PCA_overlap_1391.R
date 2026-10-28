# F1SOSS Scree Plot and 3D PCA
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
  stop("Package 'scatterplot3d' is required.")
}

# Paths
input_dir  <- "/data/02_PCA/F1SOSS/"
output_dir <- "/results/02_PCA/F1SOSS/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

eigenvec_file <- file.path(input_dir, "F1SOSS_dip_overlap_1391_PCA.eigenvec")
eigenval_file <- file.path(input_dir, "F1SOSS_dip_overlap_1391_PCA.eigenval")
pop_file      <- file.path(input_dir, "F1SOSS_pca.pop.txt")

FONT_SIZE_PT <- 14
FIGURE_DPI   <- 300

# 1. Load Data
eigvec <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
eigval <- read.table(eigenval_file, header = FALSE)
poptab <- read.table(pop_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

n_pcs <- nrow(eigval)
pcs   <- paste0("PC", seq_len(n_pcs))
percentage <- eigval$V1 / sum(eigval$V1) * 100
eigval_df  <- data.frame(PC = seq_len(n_pcs), proportion = percentage)

pca_data <- eigvec[, 2:(n_pcs + 2)]
colnames(pca_data) <- c("IID", pcs)
pca_data$IID <- as.character(pca_data$IID)

if ("vcf_id" %in% colnames(poptab)) {
  colnames(poptab)[colnames(poptab) == "vcf_id"] <- "IID"
} else if (!("IID" %in% colnames(poptab))) {
  colnames(poptab)[1] <- "IID"
}
poptab$IID <- as.character(poptab$IID)
pca_data <- merge(pca_data, poptab[, c("IID", "Group")], by = "IID", all.x = TRUE)

# ============================================================
# Figure A: Scree Plot
# ============================================================
n_scree <- min(15, n_pcs)
scree_data <- eigval_df[seq_len(n_scree), ]

p_scree <- ggplot(scree_data, aes(x = factor(PC, levels = seq_len(n_scree)), y = proportion)) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_point(size = 3, color = "#8B2671") +
  geom_line(aes(group = 1), color = "#8B2671", linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")), vjust = -0.8, size = 4) +
  labs(x = "Principal Component", y = "Variance Explained (%)") +
  theme_classic(base_size = FONT_SIZE_PT) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))

ggsave(file.path(output_dir, "F1SOSS_Scree_Plot.pdf"), p_scree, width = 8, height = 4)
ggsave(file.path(output_dir, "F1SOSS_Scree_Plot.png"), p_scree, width = 8, height = 4, dpi = FIGURE_DPI)

# ============================================================
# Figure B: 3D PCA Plot
# ============================================================
group_colors <- c("HLJ"="#8B2671", "GD"="#069F72", "SO"="#FFCA28", "SS"="#0773B2")
group_shapes <- c("HLJ"=15, "GD"=17, "SO"=16, "SS"=18)
group_order  <- c("HLJ", "GD", "SO", "SS")

save_3d_pca <- function(pca_data, percentage, dataset_name, file_stem) {
  d <- pca_data[complete.cases(pca_data[, c("PC1","PC2","PC3","Group")]), ]
  point_colors <- unname(group_colors[as.character(d$Group)])
  point_shapes <- unname(group_shapes[as.character(d$Group)])
  present_groups <- intersect(group_order, unique(as.character(d$Group)))
  
  draw_plot <- function() {
    par(mar=c(4,4,3,7), xpd=TRUE)
    scatterplot3d::scatterplot3d(
      d$PC1, d$PC2, d$PC3,
      # 【修改点】把 cex.symbols 从 1.35 调大到 2.2，点会明显变大
      color=point_colors, pch=point_shapes, cex.symbols=2.2, 
      angle=55, grid=TRUE, box=FALSE, zlim=c(-1,1),
      xlab=paste0("PC1 (", round(percentage[1],2), "%)"),
      ylab=paste0("PC2 (", round(percentage[2],2), "%)"),
      zlab=paste0("PC3 (", round(percentage[3],2), "%)"),
      main=paste0(dataset_name, " 3D PCA (PC1 + PC2 + PC3)")
    )
    legend("topright", inset=c(-0.25,0), legend=present_groups,
           col=unname(group_colors[present_groups]),
           pch=unname(group_shapes[present_groups]),
           pt.cex=1.5, # 图例里的点也稍微放大一点
           bty="n", title="Group")
  }
  
  png(paste0(file_stem, "_3D_PC1_PC2_PC3.png"), width=2400, height=2000, res=300)
  draw_plot(); dev.off()
  pdf(paste0(file_stem, "_3D_PC1_PC2_PC3.pdf"), width=8, height=7)
  draw_plot(); dev.off()
  cat("Saved: ", file_stem, "_3D_PC1_PC2_PC3.pdf/png\n", sep="")
}

# Run 3D PCA save
save_3d_pca(pca_data, percentage, "F1SOSS", file.path(output_dir, "F1SOSS_overlap_1391"))

cat("F1SOSS: Scree plot and 3D PCA saved to:", output_dir, "\n")