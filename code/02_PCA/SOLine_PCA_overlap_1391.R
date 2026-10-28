# ============================================================
# PCA plots for SOLine (overlap_1391)
# Group-coloured + Heterogeneity-coloured PC1 vs PC2 (side by side)
# + Scree plot. Font: Arial. Legend: top. NO 3D PCA.
# PDF saved with base pdf() device (Scheme A).
# Scree plot width = 2 x PCA plot width.
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ============================================================
# Paths
# ============================================================
input_dir  <- "/data/02_PCA/SOLine/"
output_dir <- "/results/02_PCA/SOLine/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================
# Output size and font settings
# ============================================================
FONT_SIZE_PT <- 13
FIGURE_DPI   <- 300

PCA_2D_WIDTH_IN  <- 6
PCA_2D_HEIGHT_IN <- 4.5

SCREE_WIDTH_IN  <- PCA_2D_WIDTH_IN * 2
SCREE_HEIGHT_IN <- PCA_2D_HEIGHT_IN

# ============================================================
# Group colours and shapes
# ============================================================
group_order <- c(
  "HLJ", "GD", "F1", "Line1", "Line2", "Line3", "Line4", "Line6",
  "Line7", "Line8", "Line10", "Line11", "Line16", "Line19", "Line24",
  "Line30", "Line31", "Line32"
)
group_colors <- setNames(c(
  "#8B2671", "#069F72", "#FFCA28", "#9BC2E6", "#F4B183", "#C15811",
  "#757575", "#F06292", "#BA124A", "#BA8CDC", "#7030A0", "#BF7F3F",
  "#0964AF", "#774F27", "#063F6E", "#68A141", "#02192C", "#385723"
), group_order)
group_shapes <- setNames(rep(16, length(group_order)), group_order)
group_shapes[c("HLJ", "GD", "F1")] <- c(15, 17, 18)

# ============================================================
# Theme
# ============================================================
base_theme <- theme_bw(base_size = FONT_SIZE_PT) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border    = element_rect(color = "black", linewidth = 0.5, fill = NA),
    text            = element_text(size = FONT_SIZE_PT),
    axis.text       = element_text(size = FONT_SIZE_PT, color = "black"),
    axis.title      = element_text(size = FONT_SIZE_PT),
    plot.title      = element_text(hjust = 0.5, size = FONT_SIZE_PT, face = "bold"),
    legend.title    = element_text(size = FONT_SIZE_PT),
    legend.text     = element_text(size = FONT_SIZE_PT),
    legend.position = "top",
    legend.box      = "horizontal",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.key.size = unit(0.5, "cm"),
    plot.margin     = unit(c(2, 2, 2, 2), "mm")
  )

# ============================================================
# Read heterogeneity data
# ============================================================
het_file_path <- file.path(input_dir, "het_1391.txt")
het_data <- read.table(het_file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("Sample", "IID", "Het_Rate") %in% colnames(het_data))) {
  stop("het_1391.txt must contain Sample, IID and Het_Rate columns")
}
het_data <- het_data[het_data$Sample == "SOLine_dip_overlap_1391", ]
if (nrow(het_data) == 0) stop("No SOLine_dip_overlap_1391 records found in het_1391.txt")
het_subset <- unique(het_data[, c("IID", "Het_Rate")])
het_subset$IID <- as.character(het_subset$IID)

# ============================================================
# Main function
# ============================================================
run_pca_plots <- function(dataset_name, eigenvec_file, eigenval_file,
                          population_file, output_prefix) {
  cat("\n========================================\nProcessing: ", dataset_name,
      "\n========================================\n", sep = "")
  
  eigvec <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
  eigval <- read.table(eigenval_file, header = FALSE)
  poptable <- read.table(population_file, header = TRUE, comment.char = "",
                         stringsAsFactors = FALSE, check.names = FALSE)
  n_pcs <- nrow(eigval)
  if (n_pcs < 2) stop(dataset_name, " has fewer than two PCs")
  
  pcs <- paste0("PC", seq_len(n_pcs))
  percentage <- eigval$V1 / sum(eigval$V1) * 100
  eigval_df <- data.frame(PC = seq_len(n_pcs), proportion = percentage)
  
  pca_data <- eigvec[, 2:(n_pcs + 2)]
  colnames(pca_data) <- c("IID", pcs)
  pca_data$IID <- as.character(pca_data$IID)
  
  if ("vcf_id" %in% colnames(poptable)) {
    colnames(poptable)[colnames(poptable) == "vcf_id"] <- "IID"
  } else if (!("IID" %in% colnames(poptable))) {
    colnames(poptable)[1] <- "IID"
  }
  if (!("Group" %in% colnames(poptable))) stop(population_file, " has no Group column")
  poptable <- unique(poptable[, c("IID", "Group")])
  poptable$IID <- as.character(poptable$IID)
  
  missing_group <- setdiff(pca_data$IID, poptable$IID)
  if (length(missing_group) > 0) warning(dataset_name, ": missing Group: ", paste(missing_group, collapse = ", "))
  pca_data <- merge(pca_data, poptable, by = "IID", all.x = TRUE)
  pca_data <- pca_data[!is.na(pca_data$Group), ]
  pca_data <- merge(pca_data, het_subset, by = "IID", all.x = TRUE)
  
  missing_het <- pca_data$IID[is.na(pca_data$Het_Rate)]
  if (length(missing_het) > 0) warning(dataset_name, ": missing Het_Rate: ", paste(missing_het, collapse = ", "))
  pca_data$Group <- factor(pca_data$Group,
                           levels = intersect(group_order, unique(pca_data$Group)))
  het_pca_data <- pca_data[!is.na(pca_data$Het_Rate), ]
  available_colors <- group_colors[names(group_colors) %in% unique(as.character(pca_data$Group))]
  available_shapes <- group_shapes[names(group_shapes) %in% unique(as.character(pca_data$Group))]
  
  axis_label <- function(pc) {
    i <- as.integer(sub("PC", "", pc))
    paste0(pc, " (", round(percentage[i], 2), "%)")
  }
  padded_range <- function(x, padding = 0.05) {
    r <- range(x, na.rm = TRUE)
    span <- diff(r)
    if (!is.finite(span) || span == 0) span <- 1
    r + c(-1, 1) * span * padding
  }
  
  limits_12 <- list(
    x = padded_range(pca_data$PC1),
    y = padded_range(pca_data$PC2)
  )
  
  # ---- Scree plot ----
  n_scree <- min(15, n_pcs)
  scree_data <- eigval_df[seq_len(n_scree), ]
  scree_plot <- ggplot(scree_data, aes(x = factor(PC, levels = seq_len(n_scree)), y = proportion)) +
    geom_col(fill = "grey50", alpha = 0.7, width = 0.6) +
    geom_point(size = 3, color = "#8B2671") +
    geom_line(aes(group = 1), color = "#8B2671", linetype = "dashed", linewidth = 0.6) +
    geom_text(aes(label = paste0(round(proportion, 1), "%")),
              vjust = -0.5, size = FONT_SIZE_PT / 2.845276) +
    labs(title = paste0(dataset_name, " PCA scree plot"),
         x = "Principal Component", y = "Variance Explained (%)") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = NULL)
  
  # ---- Group PCA ----
  p_group <- ggplot(pca_data, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
    geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.7) +
    scale_color_manual(values = available_colors, name = "Group") +
    scale_shape_manual(values = available_shapes, name = "Group") +
    labs(x = axis_label("PC1"), y = axis_label("PC2")) +
    coord_cartesian(xlim = limits_12$x, ylim = limits_12$y) +
    base_theme
  
  # ---- Heterogeneity PCA ----
  p_het <- ggplot(het_pca_data, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
    geom_point(aes(color = Het_Rate, shape = Group), size = 4, alpha = 0.7) +
    scale_color_gradient2(low = "#0773B2", mid = "#F5F5F5", high = "#8B2671",
                          midpoint = 0.5, name = "heterogeneity", limits = c(0, 1),
                          breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_shape_manual(values = available_shapes, guide = "none") +
    labs(x = axis_label("PC1"), y = axis_label("PC2")) +
    coord_cartesian(xlim = limits_12$x, ylim = limits_12$y) +
    base_theme
  
  # ============================================================
  # Save plots
  # ============================================================
  
  # ---- Scree plot ----
  ggsave(paste0(output_prefix, "_Scree_plot.png"), scree_plot,
         width = SCREE_WIDTH_IN, height = SCREE_HEIGHT_IN, units = "in",
         dpi = FIGURE_DPI, bg = "white")
  ggsave(paste0(output_prefix, "_Scree_plot.pdf"), scree_plot,
         width = SCREE_WIDTH_IN, height = SCREE_HEIGHT_IN, units = "in",
         device = pdf, bg = "white")
  cat("Saved: ", output_prefix, "_Scree_plot.pdf/png\n", sep = "")
  
  # ---- Group PCA ----
  ggsave(paste0(output_prefix, "_PC1_PC2_Group.png"), p_group,
         width = PCA_2D_WIDTH_IN, height = PCA_2D_HEIGHT_IN, units = "in",
         dpi = FIGURE_DPI, bg = "white")
  ggsave(paste0(output_prefix, "_PC1_PC2_Group.pdf"), p_group,
         width = PCA_2D_WIDTH_IN, height = PCA_2D_HEIGHT_IN, units = "in",
         device = pdf, bg = "white")
  
  # ---- Heterogeneity PCA ----
  ggsave(paste0(output_prefix, "_PC1_PC2_Heterogeneity.png"), p_het,
         width = PCA_2D_WIDTH_IN, height = PCA_2D_HEIGHT_IN, units = "in",
         dpi = FIGURE_DPI, bg = "white")
  ggsave(paste0(output_prefix, "_PC1_PC2_Heterogeneity.pdf"), p_het,
         width = PCA_2D_WIDTH_IN, height = PCA_2D_HEIGHT_IN, units = "in",
         device = pdf, bg = "white")
  
  # ---- Combined side-by-side ----
  combined <- p_group + p_het + plot_layout(ncol = 2)
  ggsave(paste0(output_prefix, "_PC1_PC2_Group_Heterogeneity_combined.png"), combined,
         width = PCA_2D_WIDTH_IN * 1.5, height = PCA_2D_HEIGHT_IN, units = "in",
         dpi = FIGURE_DPI, bg = "white")
  ggsave(paste0(output_prefix, "_PC1_PC2_Group_Heterogeneity_combined.pdf"), combined,
         width = PCA_2D_WIDTH_IN * 1.5, height = PCA_2D_HEIGHT_IN, units = "in",
         device = pdf, bg = "white")
  
  cat("Saved: individual and combined PCA plots for ", dataset_name, "\n", sep = "")
  cat("Total PCs:", n_pcs, "\nSamples retained:", nrow(pca_data),
      "\nSamples with Heterogeneity:", nrow(het_pca_data), "\n")
  print(eigval_df[seq_len(min(5, n_pcs)), ])
}

# ============================================================
# Run
# ============================================================
run_pca_plots(
  dataset_name    = "SOLine",
  eigenvec_file   = file.path(input_dir, "SOLine_dip_overlap_1391_PCA.eigenvec"),
  eigenval_file   = file.path(input_dir, "SOLine_dip_overlap_1391_PCA.eigenval"),
  population_file = file.path(input_dir, "SOLine_pca.pop.txt"),
  output_prefix   = file.path(output_dir, "SOLine_overlap_1391")
)

cat("\nSOLine PCA plots completed.\n")