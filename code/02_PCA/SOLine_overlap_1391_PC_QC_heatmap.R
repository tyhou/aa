# PCA vs QC/covariate association heatmap
suppressPackageStartupMessages({
  library(readxl)
  library(pheatmap)
})

# Paths
input_dir  <- "/data/02_PCA/SOLine/"
output_dir <- "/results/02_PCA/SOLine/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

eigenvec_file <- file.path(input_dir, "SOLine_dip_overlap_1391_PCA.eigenvec")
eigenval_file <- file.path(input_dir, "SOLine_dip_overlap_1391_PCA.eigenval")
metadata_file <- file.path("/data/02_PCA/PCA_QC_metadata.xlsx")

# Settings
FONT_FAMILY    <- "Arial"
FONT_SIZE_PT   <- 14
FIGURE_DPI     <- 300
N_PCS_TO_CHECK <- 10
heat_vars      <- c("Depth", "MappingRate", "Heterogeneity", "Group", "DPI", "Batch")

# Read data
eigvec <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
eigval <- read.table(eigenval_file, header = FALSE)
meta   <- as.data.frame(read_excel(metadata_file, sheet = 1))

n_pcs <- nrow(eigval)
pcs   <- paste0("PC", seq_len(n_pcs))

pca_data <- eigvec[, 2:(n_pcs + 2)]
colnames(pca_data) <- c("IID", pcs)
pca_data$IID <- as.character(pca_data$IID)
meta$IID     <- as.character(meta$IID)

if (length(setdiff(heat_vars, colnames(meta))) > 0) {
  stop("Metadata missing columns: ", paste(setdiff(heat_vars, colnames(meta)), collapse = ", "))
}

merged <- merge(pca_data, meta, by = "IID", all = FALSE)
cat("Samples in PCA:", nrow(pca_data), "| in meta:", nrow(meta), "| merged:", nrow(merged), "\n")

# Encode covariates
parse_dpi <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(sub("^[Dd]", "", trimws(as.character(x)))))
}

merged$DPI_num   <- parse_dpi(merged$DPI)
merged$Batch_num <- if (is.numeric(merged$Batch)) as.numeric(merged$Batch) else as.numeric(factor(merged$Batch))
merged$Group_num <- as.numeric(factor(merged$Group))

covar_map <- c(Depth = "Depth", MappingRate = "MappingRate", Heterogeneity = "Heterogeneity",
               Group = "Group_num", DPI = "DPI_num", Batch = "Batch_num")

# Calculate correlations
pcs_check <- pcs[seq_len(min(N_PCS_TO_CHECK, n_pcs))]
rho_mat   <- matrix(NA_real_, length(pcs_check), length(heat_vars), dimnames = list(pcs_check, heat_vars))
pval_mat  <- rho_mat

for (pc in pcs_check) {
  for (v in heat_vars) {
    x  <- merged[[pc]]
    y  <- merged[[covar_map[[v]]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 5) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      rho_mat[pc, v]  <- unname(ct$estimate)
      pval_mat[pc, v] <- ct$p.value
    }
  }
}

# Save table
out_df <- data.frame(
  PC          = rep(rownames(rho_mat), times = ncol(rho_mat)),
  Variable    = rep(colnames(rho_mat), each  = nrow(rho_mat)),
  SpearmanRho = as.vector(rho_mat),
  P_value     = as.vector(pval_mat),
  stringsAsFactors = FALSE
)
write.csv(out_df, file.path(output_dir, "SOLine_overlap_1391_PC_QC_heatmap.csv"), row.names = FALSE)

# Generate heatmap
star <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
disp <- matrix(sprintf("%.2f%s", rho_mat, star(pval_mat)), nrow = nrow(rho_mat), dimnames = dimnames(rho_mat))

heat_palette <- colorRampPalette(c("#0773B2", "white", "#8B2671"))(101)
breaks_seq   <- seq(-1, 1, length.out = 102)

draw_heatmap <- function() {
  pheatmap(
    rho_mat, display_numbers = disp, number_color = "black",
    fontsize = FONT_SIZE_PT, fontsize_number = FONT_SIZE_PT - 2,
    fontfamily = FONT_FAMILY, color = heat_palette, breaks = breaks_seq,
    cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey80", angle_col = 45,
    main = "PC association with QC metrics and batch variables"
  )
}

png(file.path(output_dir, "SOLine_overlap_1391_PC_QC_heatmap.png"), width = 7, height = 5.5, units = "in", res = FIGURE_DPI)
draw_heatmap(); dev.off()

cairo_pdf(file.path(output_dir, "SOLine_overlap_1391_PC_QC_heatmap.pdf"), width = 7, height = 5.5, family = FONT_FAMILY)
draw_heatmap(); dev.off()

cat("Heatmap and CSV saved to:", output_dir, "\n")