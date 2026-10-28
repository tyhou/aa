# ============================================================
# F1SOSS Allele Frequency Maps
# ============================================================

library(ggplot2)

# ---- Paths ----
input_dir  <- "/data/04_Allele_Frequency/F1SOSS"
output_dir <- "/results/04_Allele_Frequency/F1SOSS"

# Shared helper functions live in this directory; run.sh exports the path
code_dir   <- Sys.getenv("ALLELEFREQ_CODE_DIR", unset = "/code/04_Allele_Frequency")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Source shared functions
source(file.path(code_dir, "allelefreq_functions_1391.R"))

# ---- Read VCF ----
vcf_file <- file.path(input_dir, "F1SOSS_hap_overlap_1391.vcf")
vcf <- read_haploid_vcf(vcf_file)

# ---- Identify informative sites ----
mask <- informative_mask(vcf, "11730", "22971")
cat("F1SOSS VCF records:", nrow(vcf$fixed), "\nInformative sites:", sum(mask, na.rm=TRUE), "\n")

# ---- Make frequency table ----
freq <- make_frequency_table(vcf, mask, "F1SOSS", out_dir = output_dir)

# ---- Plot maps ----
f1_groups <- list(
  F1_SO = c("W585", "W586", "W587", "W588", "W589", "W590",
            "W591", "W592", "W594", "W595", "W596"),
  F1_SS = c("W680", "W681", "W682", "W788", "W789", "W790",
            "W791", "W792", "W793", "W794", "W795"),
  Selected_W586_W681_W682 = c("W586", "W681", "W682"),
  Remaining_15_samples = c(
    "W585", "W587", "W588", "W589", "W590",
    "W591", "W592", "W594", "W595", "W596",
    "W680", "W791", "W792", "W793", "W795"
  )
)

plot_frequency_maps(freq, vcf$contigs, "F1SOSS", "F1SOSS", 
                    combined_groups = f1_groups, 
                    out_dir = output_dir)

cat("F1SOSS allele-frequency maps completed.\n")