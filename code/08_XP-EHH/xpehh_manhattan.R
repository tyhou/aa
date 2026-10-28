library(rehh)
library(dplyr)
library(ggplot2)
library(stringr)

# Paths and settings
input_dir  <- "/data/08_XP-EHH"
output_dir <- "/results/08_XP-EHH"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gff_file <- file.path(input_dir, "IId_11730_genome.gff3")
chromosomes <- 1:8
p_threshold <- 0.05
merge_distance_bp <- 200000
font_size <- 14

# Read VCF and calculate iES statistics
read_and_scan <- function(group, chr) {
  input_file <- file.path(input_dir, "phasing", sprintf("phasing_%s_chr%d.vcf.gz", group, chr))
  hap <- data2haplohh(
    hap_file = input_file,
    polarize_vcf = FALSE,
    vcf_reader = "data.table",
    remove_multiple_markers = TRUE
  )
  scan_hh(hap, discard_integration_at_border = FALSE)
}

scan_by_group <- list(high = NULL, low = NULL, all = NULL)

for (chr in chromosomes) {
  for (group in names(scan_by_group)) {
    chromosome_scan <- read_and_scan(group, chr)
    scan_by_group[[group]] <- bind_rows(scan_by_group[[group]], chromosome_scan)
  }
  cat(sprintf(
    "Completed Chr%d: high=%d, low=%d, all=%d markers\n",
    chr,
    sum(scan_by_group$high$CHR == chr),
    sum(scan_by_group$low$CHR == chr),
    sum(scan_by_group$all$CHR == chr)
  ))
}

check_chromosomes <- function(scan_data, group_name) {
  observed <- sort(unique(as.integer(as.character(scan_data$CHR))))
  if (!identical(observed, chromosomes)) {
    stop(group_name, " is incomplete. Observed chromosomes: ", paste(observed, collapse = ", "))
  }
}

invisible(Map(check_chromosomes, scan_by_group, names(scan_by_group)))

# Calculate genome-wide iHS and XP-EHH
ihs_high <- ihh2ihs(scan_by_group$high)[["ihs"]]
ihs_low <- ihh2ihs(scan_by_group$low)[["ihs"]]
ihs_all <- ihh2ihs(scan_by_group$all)[["ihs"]]

xpehh_high_low <- ies2xpehh(scan_by_group$high, scan_by_group$low, "high", "low")
xpehh_high_all <- ies2xpehh(scan_by_group$high, scan_by_group$all, "high", "all")

write.table(ihs_high, file.path(output_dir, "ihs_high.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ihs_low, file.path(output_dir, "ihs_low.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ihs_all, file.path(output_dir, "ihs_all.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(xpehh_high_low, file.path(output_dir, "xpehh_high_vs_low.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(xpehh_high_all, file.path(output_dir, "xpehh_high_vs_all.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Annotate significant XP-EHH sites and selection direction
logp_threshold <- -log10(p_threshold)
xpehh_annotated <- xpehh_high_low %>%
  mutate(
    P_VALUE = 10^(-LOGPVALUE),
    IS_SIGNIFICANT = LOGPVALUE > logp_threshold,
    SELECTION_DIRECTION = case_when(
      IS_SIGNIFICANT & XPEHH_high_low > 0 ~ "high",
      IS_SIGNIFICANT & XPEHH_high_low < 0 ~ "low",
      TRUE ~ "not_significant"
    ),
    SNP_ID = paste0("Chr", CHR, "_", POSITION)
  )

write.table(
  xpehh_annotated, file.path(output_dir, "xpehh_high_vs_low_annotated.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Merge nearby significant sites into candidate regions
merge_candidate_regions <- function(data, distance_bp) {
  sites <- data %>%
    filter(IS_SIGNIFICANT) %>%
    arrange(as.integer(CHR), POSITION)
  
  if (nrow(sites) == 0) return(data.frame())
  
  region_number <- integer(nrow(sites))
  region_number[1] <- 1L
  for (i in 2:nrow(sites)) {
    new_region <- sites$CHR[i] != sites$CHR[i - 1] ||
      sites$POSITION[i] - sites$POSITION[i - 1] > distance_bp
    region_number[i] <- region_number[i - 1] + as.integer(new_region)
  }
  
  sites %>%
    mutate(REGION_NUMBER = region_number) %>%
    group_by(REGION_NUMBER, CHR) %>%
    summarise(
      START = min(POSITION),
      END = max(POSITION),
      REGION_LENGTH_BP = END - START + 1,
      N_SITES = n(),
      DIRECTION = if_else(
        n_distinct(SELECTION_DIRECTION) == 1,
        first(SELECTION_DIRECTION),
        "mixed"
      ),
      MAX_LOGPVALUE = max(LOGPVALUE),
      MEAN_XPEHH = mean(XPEHH_high_low),
      LEAD_SNP = SNP_ID[which.max(LOGPVALUE)],
      LEAD_POSITION = POSITION[which.max(LOGPVALUE)],
      SITE_LIST = paste(SNP_ID, collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      REGION_ID = sprintf("Region_%02d_Chr%s", row_number(), CHR),
      .before = 1
    ) %>%
    select(-REGION_NUMBER)
}

candidate_regions <- merge_candidate_regions(xpehh_annotated, merge_distance_bp)
write.table(
  candidate_regions, file.path(output_dir, "xpehh_candidate_regions.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Read GFF3 and build gene annotations
extract_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  value <- str_match(attributes, pattern)[, 2]
  ifelse(is.na(value), NA_character_, URLdecode(value))
}

gff <- read.delim(
  gff_file,
  sep = "\t", header = FALSE, comment.char = "#",
  quote = "", stringsAsFactors = FALSE,
  col.names = c("CHR", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRIBUTES")
) %>%
  mutate(
    CHR = as.character(CHR),
    FEATURE_ID = extract_attribute(ATTRIBUTES, "ID"),
    LOCUS_TAG = extract_attribute(ATTRIBUTES, "locus_tag"),
    GENE_NAME = extract_attribute(ATTRIBUTES, "gene"),
    GENE_BIOTYPE = extract_attribute(ATTRIBUTES, "gene_biotype"),
    PRODUCT_LOWER = extract_attribute(ATTRIBUTES, "product"),
    PRODUCT_UPPER = extract_attribute(ATTRIBUTES, "Product"),
    PRODUCT = coalesce(PRODUCT_LOWER, PRODUCT_UPPER)
  )

gene_products <- gff %>%
  filter(!is.na(LOCUS_TAG), !is.na(PRODUCT)) %>%
  group_by(LOCUS_TAG) %>%
  summarise(PRODUCT = first(PRODUCT), .groups = "drop")

genes <- gff %>%
  filter(TYPE == "gene") %>%
  transmute(
    CHR,
    GENE_START = START,
    GENE_END = END,
    STRAND,
    GENE_ID = coalesce(LOCUS_TAG, str_remove(FEATURE_ID, "^gene-")),
    GENE_NAME,
    GENE_BIOTYPE
  ) %>%
  left_join(gene_products, by = c("GENE_ID" = "LOCUS_TAG")) %>%
  arrange(as.integer(CHR), GENE_START)

# Map positive significant sites to directly overlapping and nearest genes
high_candidate_sites <- xpehh_annotated %>%
  filter(IS_SIGNIFICANT, XPEHH_high_low > 0) %>%
  arrange(as.integer(CHR), POSITION)

annotate_site <- function(site_row, gene_table) {
  chr_genes <- gene_table %>% filter(CHR == as.character(site_row$CHR))
  overlapping <- chr_genes %>%
    filter(GENE_START <= site_row$POSITION, GENE_END >= site_row$POSITION)
  
  common <- data.frame(
    CHR = site_row$CHR,
    POSITION = site_row$POSITION,
    SNP_ID = site_row$SNP_ID,
    XPEHH_high_low = site_row$XPEHH_high_low,
    LOGPVALUE = site_row$LOGPVALUE,
    P_VALUE = site_row$P_VALUE,
    stringsAsFactors = FALSE
  )
  
  if (nrow(overlapping) > 0) {
    return(bind_cols(
      common[rep(1, nrow(overlapping)), ],
      data.frame(
        ANNOTATION = "gene_overlap",
        GENE_ID = overlapping$GENE_ID,
        GENE_NAME = overlapping$GENE_NAME,
        GENE_START = overlapping$GENE_START,
        GENE_END = overlapping$GENE_END,
        STRAND = overlapping$STRAND,
        GENE_BIOTYPE = overlapping$GENE_BIOTYPE,
        PRODUCT = overlapping$PRODUCT,
        NEAREST_GENE_ID = NA_character_,
        DISTANCE_TO_NEAREST_GENE_BP = 0,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  distance <- pmax(chr_genes$GENE_START - site_row$POSITION, 0) +
    pmax(site_row$POSITION - chr_genes$GENE_END, 0)
  nearest_index <- which.min(distance)
  nearest <- chr_genes[nearest_index, ]
  
  bind_cols(
    common,
    data.frame(
      ANNOTATION = "intergenic",
      GENE_ID = NA_character_,
      GENE_NAME = NA_character_,
      GENE_START = NA_integer_,
      GENE_END = NA_integer_,
      STRAND = NA_character_,
      GENE_BIOTYPE = NA_character_,
      PRODUCT = NA_character_,
      NEAREST_GENE_ID = nearest$GENE_ID,
      DISTANCE_TO_NEAREST_GENE_BP = distance[nearest_index],
      stringsAsFactors = FALSE
    )
  )
}

site_gene_map <- if (nrow(high_candidate_sites) > 0) {
  bind_rows(lapply(seq_len(nrow(high_candidate_sites)), function(i) {
    annotate_site(high_candidate_sites[i, ], genes)
  }))
} else {
  data.frame()
}

write.table(
  site_gene_map, file.path(output_dir, "xpehh_high_candidate_site_gene_map.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# Select genes by candidate-region membership
annotate_region_snp_genes <- function(region_row, gene_table, snp_table) {
  region_genes <- gene_table %>%
    filter(
      CHR == as.character(region_row$CHR),
      GENE_START <= region_row$END,
      GENE_END >= region_row$START
    )
  
  if (nrow(region_genes) == 0) return(data.frame())
  
  annotated <- lapply(seq_len(nrow(region_genes)), function(i) {
    gene <- region_genes[i, ]
    gene_snps <- snp_table %>%
      filter(
        CHR == as.character(region_row$CHR),
        POSITION >= max(gene$GENE_START, region_row$START),
        POSITION <= min(gene$GENE_END, region_row$END)
      ) %>%
      arrange(POSITION)
    
    if (nrow(gene_snps) == 0) return(NULL)
    
    lead_index <- which.max(gene_snps$LOGPVALUE)
    data.frame(
      REGION_ID = region_row$REGION_ID,
      CHR = as.character(region_row$CHR),
      REGION_START = region_row$START,
      REGION_END = region_row$END,
      REGION_LENGTH_BP = region_row$REGION_LENGTH_BP,
      REGION_SIGNIFICANT_SITE_COUNT = region_row$N_SITES,
      REGION_DIRECTION = region_row$DIRECTION,
      REGION_MAX_LOGPVALUE = region_row$MAX_LOGPVALUE,
      REGION_LEAD_SNP = region_row$LEAD_SNP,
      GENE_ID = gene$GENE_ID,
      GENE_NAME = gene$GENE_NAME,
      GENE_START = gene$GENE_START,
      GENE_END = gene$GENE_END,
      STRAND = gene$STRAND,
      GENE_BIOTYPE = gene$GENE_BIOTYPE,
      PRODUCT = gene$PRODUCT,
      SNP_COUNT = nrow(gene_snps),
      SIGNIFICANT_SNP_COUNT = sum(gene_snps$IS_SIGNIFICANT),
      HIGH_SIGNIFICANT_SNP_COUNT = sum(
        gene_snps$IS_SIGNIFICANT & gene_snps$XPEHH_high_low > 0
      ),
      LOW_SIGNIFICANT_SNP_COUNT = sum(
        gene_snps$IS_SIGNIFICANT & gene_snps$XPEHH_high_low < 0
      ),
      FIRST_SNP_POSITION = min(gene_snps$POSITION),
      LAST_SNP_POSITION = max(gene_snps$POSITION),
      MAX_GENE_LOGPVALUE = gene_snps$LOGPVALUE[lead_index],
      GENE_LEAD_SNP = gene_snps$SNP_ID[lead_index],
      GENE_LEAD_XPEHH = gene_snps$XPEHH_high_low[lead_index],
      SNP_LIST = paste(gene_snps$SNP_ID, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(annotated)
}

region_snp_genes <- if (nrow(candidate_regions) > 0) {
  bind_rows(lapply(seq_len(nrow(candidate_regions)), function(i) {
    annotate_region_snp_genes(candidate_regions[i, ], genes, xpehh_annotated)
  })) %>%
    arrange(as.integer(CHR), REGION_START, GENE_START)
} else {
  data.frame()
}

write.table(
  region_snp_genes, file.path(output_dir, "xpehh_candidate_region_snp_genes.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# High-virulence gene summary
gene_summary <- region_snp_genes %>%
  filter(REGION_DIRECTION == "high")

write.table(
  gene_summary, file.path(output_dir, "xpehh_high_candidate_gene_summary.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# Build chromosome length table
chromosome_layout <- gff %>%
  filter(CHR %in% as.character(chromosomes)) %>%
  group_by(CHR) %>%
  summarise(CHR_LENGTH = max(END), .groups = "drop") %>%
  mutate(
    CHR_ORDER = as.integer(CHR),
    CHR = factor(paste0("Chr", CHR), levels = paste0("Chr", chromosomes))
  ) %>%
  arrange(CHR_ORDER)

# Plot data
plot_data <- xpehh_annotated %>%
  mutate(
    CHR = factor(paste0("Chr", CHR), levels = paste0("Chr", chromosomes)),
    SIGNED_LOGPVALUE = if_else(XPEHH_high_low >= 0, LOGPVALUE, -LOGPVALUE)
  )

# Dummy points to force x-axis span (Fixed for dplyr >= 1.1.0)
dummy_points <- chromosome_layout %>%
  group_by(CHR) %>%
  reframe(
    POSITION = c(0, CHR_LENGTH),
    SIGNED_LOGPVALUE = c(0, 0),
    IS_SIGNIFICANT = FALSE,
    SELECTION_DIRECTION = NA_character_
  )

plot_data_padded <- plot_data %>%
  select(CHR, POSITION, SIGNED_LOGPVALUE, IS_SIGNIFICANT, SELECTION_DIRECTION) %>%
  bind_rows(dummy_points) %>%
  mutate(CHR = factor(CHR, levels = levels(plot_data$CHR)))

y_limit <- max(abs(plot_data$SIGNED_LOGPVALUE), na.rm = TRUE)
y_limit <- ceiling((y_limit * 1.08) * 2) / 2

direction_colors <- c(high = "#D73027", low = "#225EA8")

# Point aesthetics
point_size  <- 4.0
point_alpha <- 0.5

xpehh_plot <- ggplot(plot_data_padded, aes(POSITION, SIGNED_LOGPVALUE)) +
  geom_blank(aes(x = POSITION, y = 0)) +
  geom_point(
    data = plot_data,
    color = "gray30", shape = 16,
    size = point_size, alpha = point_alpha
  ) +
  geom_point(
    data = filter(plot_data, IS_SIGNIFICANT),
    aes(fill = SELECTION_DIRECTION),
    shape = 21, color = "white", stroke = 0.3,
    size = point_size, alpha = point_alpha
  ) +
  geom_hline(
    yintercept = c(-logp_threshold, logp_threshold),
    linewidth = 0.55, linetype = "dashed", color = "#444444"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#666666") +
  scale_fill_manual(
    values = direction_colors,
    breaks = c("high", "low"),
    labels = c("High virulence", "Low virulence"),
    name = "Selection direction"
  ) +
  scale_x_continuous(
    labels = function(x) format(x / 1e6, nsmall = 1, trim = TRUE),
    expand = expansion(mult = c(0.02, 0.02)),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_y_continuous(
    breaks = pretty(c(-y_limit, y_limit), n = 8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(ylim = c(-y_limit, y_limit), clip = "off") +
  facet_wrap(~ CHR, nrow = 1, scales = "free_x") +
  labs(
    title = "XP-EHHL high vs low (Fws > 0.5)",
    x = "Position (Mb)",
    y = expression(paste("Signed -log"[10], "(", italic(P), ")"))
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "white",
                          stroke = 0.3, size = 3)
    )
  ) +
  theme_bw(base_size = font_size) +
  theme(
    text = element_text(size = font_size, color = "black"),
    plot.title = element_text(size = font_size, face = "bold", hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size, color = "black"),
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.position = c(0.995,0.1),
    legend.justification = c(1, 0.5),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), color = NA),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(0.15, "cm"),
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.5),
    strip.text = element_text(size = font_size, face = "bold"),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(output_dir, "xpehh_high_vs_low_manhattan.pdf"),
  plot = xpehh_plot, width = 15, height = 4, units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "xpehh_high_vs_low_manhattan.png"),
  plot = xpehh_plot, width = 15, height = 4, units = "in",
  dpi = 300, bg = "white"
)

# Write report
report_lines <- c(
  "XP-EHH high-versus-low analysis",
  sprintf("Chromosomes analyzed: %s", paste(chromosomes, collapse = ", ")),
  sprintf("XP-EHH sites: %d", nrow(xpehh_annotated)),
  sprintf("Significance threshold: P < %.3f (-log10(P) > %.3f)", p_threshold, logp_threshold),
  sprintf("Significant sites: %d", sum(xpehh_annotated$IS_SIGNIFICANT)),
  sprintf("High-virulence candidate sites: %d", nrow(high_candidate_sites)),
  sprintf("Low-virulence candidate sites: %d", sum(xpehh_annotated$IS_SIGNIFICANT & xpehh_annotated$XPEHH_high_low < 0)),
  sprintf("Candidate regions: %d", nrow(candidate_regions)),
  sprintf("Candidate-region genes containing analyzed SNPs: %d", nrow(region_snp_genes)),
  sprintf("High-virulence-region genes containing analyzed SNPs: %d", nrow(gene_summary)),
  sprintf("Intergenic high-virulence candidate sites: %d", sum(site_gene_map$ANNOTATION == "intergenic")),
  "Positive XP-EHH indicates longer haplotype homozygosity in the high-virulence group relative to the low-virulence group."
)
writeLines(report_lines, file.path(output_dir, "xpehh_analysis_report.txt"))

# Save workspace
save(
  scan_by_group, ihs_high, ihs_low, ihs_all,
  xpehh_high_low, xpehh_high_all, xpehh_annotated,
  candidate_regions, genes, high_candidate_sites,
  site_gene_map, region_snp_genes, gene_summary, chromosome_layout,
  file = file.path(output_dir, "xpehh_analysis.RData")
)

cat("\nAnalysis completed.\n")
cat("XP-EHH sites:", nrow(xpehh_annotated), "\n")
cat("Significant sites:", sum(xpehh_annotated$IS_SIGNIFICANT), "\n")
cat("High-virulence candidate sites:", nrow(high_candidate_sites), "\n")
cat("Candidate-region genes containing analyzed SNPs:", nrow(region_snp_genes), "\n")
cat("High-virulence-region genes containing analyzed SNPs:", nrow(gene_summary), "\n")