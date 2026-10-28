# Load packages
library(tidyverse)
library(writexl)
library(ggplot2)
library(gridExtra)

# Paths
INPUT_DIR  <- "/data/06_Fst_Pi"
OUTPUT_DIR <- "/results/06_Fst_Pi"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

OUTPUT_XLSX <- file.path(OUTPUT_DIR, "merge_pi_Fst-1000win-100step.xlsx")
GFF_FILE    <- file.path(INPUT_DIR, "IId_11730_genome.gff3")

# Figure settings (Removed specific font family to use default sans)
FONT_SIZE_PT <- 18
FIGURE_WIDTH_IN <- 8
FIGURE_HEIGHT_IN <- 8
MAIN_PANEL_ASPECT_RATIO <- 1
FIGURE_DPI <- 300

# Read data
dt_low <- read.delim(file.path(INPUT_DIR, "SOLine_dip_overlap_1391_low_pi_1000win_100step.windowed.pi"), sep="\t", header=TRUE, check.names=FALSE)
dt_high <- read.delim(file.path(INPUT_DIR, "SOLine_dip_overlap_1391_high_pi_1000win_100step.windowed.pi"), sep="\t", header=TRUE, check.names=FALSE)
fst_high_low <- read.delim(file.path(INPUT_DIR, "SOLine_dip_overlap_1391_fst_1000win_100step.windowed.weir.fst"), sep="\t", header=TRUE, check.names=FALSE)

# Merge data
dt2 <- merge(dt_low, dt_high, by = c("CHROM","BIN_START","BIN_END"), all = TRUE)
dt3 <- merge(dt2, fst_high_low, by = c("CHROM","BIN_START","BIN_END"), all = TRUE)

names(dt3)[names(dt3) =="PI.x"] <- "PI.low"
names(dt3)[names(dt3) =="PI.y"] <- "PI.high"
dt3[is.na(dt3)] <- 0

df3 <- dt3[,c("CHROM","BIN_START","BIN_END","PI.low","PI.high","N_VARIANTS","WEIGHTED_FST")]
df3$ROD <- df3$PI.high / df3$PI.low
df3$log2ROD <- log2(df3$ROD)
df3 <- df3[!is.infinite(df3$log2ROD),]

# Calculate thresholds
vline1 <- quantile(df3$log2ROD, 0.75)
vline2 <- quantile(df3$log2ROD, 0.25)
hline <- quantile(df3$WEIGHTED_FST, 0.75)

# Select regions
select_gene_low <- df3[df3$log2ROD > vline1 & df3$WEIGHTED_FST > hline, ]
select_gene_high <- df3[df3$log2ROD < vline2 & df3$WEIGHTED_FST > hline, ]

# Parse GFF3
extract_gff_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]*)")
  value <- stringr::str_match(attributes, pattern)[, 2]
  value[value == ""] <- NA_character_
  vapply(value, function(x) {
    if (is.na(x)) NA_character_ else utils::URLdecode(x)
  }, character(1))
}

gff <- read.delim(GFF_FILE, sep = "\t", header = FALSE, comment.char = "#", quote = "", fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
gff <- gff[, seq_len(9)]
names(gff) <- c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")

gff$locus_tag <- extract_gff_attribute(gff$attributes, "locus_tag")
gff$Name <- extract_gff_attribute(gff$attributes, "Name")
gff$gene <- extract_gff_attribute(gff$attributes, "gene")
gff$Product <- extract_gff_attribute(gff$attributes, "Product")
gff$product <- extract_gff_attribute(gff$attributes, "product")
gff$record_name <- dplyr::coalesce(gff$Name, gff$gene, gff$Product, gff$product)

gene_names <- gff |>
  dplyr::filter(!is.na(locus_tag), feature %in% c("mRNA", "rRNA", "tRNA", "ncRNA"), !is.na(record_name)) |>
  dplyr::group_by(locus_tag) |>
  dplyr::summarise(gene_name = dplyr::first(record_name), .groups = "drop")

genes <- gff |>
  dplyr::filter(feature == "gene", !is.na(locus_tag)) |>
  dplyr::transmute(seqid = as.character(seqid), gene_start = as.numeric(start), gene_end = as.numeric(end), locus_tag, gene_row_name = dplyr::coalesce(Name, gene)) |>
  dplyr::left_join(gene_names, by = "locus_tag") |>
  dplyr::mutate(gene_name = dplyr::coalesce(gene_row_name, gene_name, locus_tag)) |>
  dplyr::arrange(seqid, gene_start, gene_end)

annotate_gene_overlaps <- function(windows, gene_table) {
  windows$gene_locus_tag <- NA_character_
  windows$gene_name <- NA_character_
  for (i in seq_len(nrow(windows))) {
    hits <- gene_table[
      gene_table$seqid == as.character(windows$CHROM[i]) &
      gene_table$gene_start <= windows$BIN_END[i] &
      gene_table$gene_end >= windows$BIN_START[i],
    ]
    if (nrow(hits) > 0) {
      windows$gene_locus_tag[i] <- paste(hits$locus_tag, collapse = "; ")
      windows$gene_name[i] <- paste(hits$gene_name, collapse = "; ")
    }
  }
  windows
}

select_gene_high <- annotate_gene_overlaps(select_gene_high, genes)

# Save to Excel
xl_list <- list()
xl_list$merge_pi_Fst <- dplyr::select(df3, -BIN_END)
xl_list$select_gene_low <- dplyr::select(select_gene_low, -BIN_END)
xl_list$select_gene_high <- dplyr::select(select_gene_high, -BIN_END)
write_xlsx(xl_list, OUTPUT_XLSX, format_headers = FALSE)

# Plot main scatter
p <- ggplot(df3, aes(log2ROD, WEIGHTED_FST)) +
  annotate("rect", xmin = -Inf, xmax = vline2, ymin = hline, ymax = Inf, fill = "red", alpha = 0.16) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = hline, color = "red", linewidth = 1, linetype = "dotted") + 
  geom_vline(xintercept = vline1, color = "red", linewidth = 1, linetype = "dotted") + 
  geom_vline(xintercept = vline2, color = "red", linewidth = 1, linetype = "dotted") + 
  scale_color_brewer(palette = "Set1") +
  labs(x = expression(log[2]("Pi"[High]/"Pi"[Low])), y = "Fst")

mytheme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  text = element_text(size = FONT_SIZE_PT),
  legend.title = element_blank(),
  legend.position = c(0.9, 1.02),
  legend.justification = c("top"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.box.background = element_rect(fill = NA, colour = NA),
  axis.line.y = element_line(color = "black"),
  axis.line.x = element_line(color = "black"),
  axis.title = element_text(size = FONT_SIZE_PT, color = "black"),
  axis.text = element_text(size = FONT_SIZE_PT, color = "black"),
  aspect.ratio = MAIN_PANEL_ASPECT_RATIO,
  plot.margin = unit(x = c(0, 0, 0, 0), units = "inches")
)
p <- p + theme_bw(base_size = FONT_SIZE_PT) + mytheme

# Plot marginal histograms
p1 <- ggplot(df3) +
  geom_histogram(aes(x = log2ROD, y = after_stat(density)), fill = '#dedede', colour = "#dedede", binwidth = 0.08) +
  theme_classic(base_size = FONT_SIZE_PT) +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(df3) +
  geom_histogram(aes(x = after_stat(density), y = WEIGHTED_FST), fill = '#dedede', colour = "#dedede", binwidth = 0.01) +
  theme_classic(base_size = FONT_SIZE_PT) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

empty <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p3 <- grid.arrange(p1, empty, p, p2, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

# Save plots (Removed family argument to avoid font errors)
ggsave(file.path(OUTPUT_DIR, "ROD_Fst_low_high_1000win_100step.pdf"), plot = p, device = grDevices::cairo_pdf, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in", dpi = FIGURE_DPI, bg = "white")
ggsave(file.path(OUTPUT_DIR, "ROD_Fst_low_high_1000win_100step.png"), plot = p, device = "png", width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in", dpi = FIGURE_DPI, bg = "white")
ggsave(file.path(OUTPUT_DIR, "ROD_Fst_low_high_1000win_100step_density.pdf"), plot = p3, device = grDevices::cairo_pdf, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in", dpi = FIGURE_DPI, bg = "white")
ggsave(file.path(OUTPUT_DIR, "ROD_Fst_low_high_1000win_100step_density.png"), plot = p3, device = "png", width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in", dpi = FIGURE_DPI, bg = "white")