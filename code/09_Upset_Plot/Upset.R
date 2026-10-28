# Load libraries
library(UpSetR)
library(ggplot2)
library(openxlsx)
library(dplyr)

# Paths
input_dir  <- "/data/09_Upset_Plot"
output_dir <- "/results/09_Upset_Plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read gene lists
list1 <- scan(file.path(input_dir, "list1-3-8-11-24-genes.txt"), what = character(), quiet = TRUE)
list2 <- scan(file.path(input_dir, "list2-xpehh-genes.txt"), what = character(), quiet = TRUE)
list3 <- scan(file.path(input_dir, "list3-fst-genes.txt"), what = character(), quiet = TRUE)
list4 <- scan(file.path(input_dir, "list4-highlyPolymorphic-genes.txt"), what = character(), quiet = TRUE)

# Build binary matrix
all_unique_genes <- unique(c(list1, list2, list3, list4))

binary_matrix <- data.frame(
  Gene = all_unique_genes,
  List1 = as.integer(all_unique_genes %in% list1),
  List2 = as.integer(all_unique_genes %in% list2),
  List3 = as.integer(all_unique_genes %in% list3),
  List4 = as.integer(all_unique_genes %in% list4)
)

# Prepare data for UpSet plot
upset_data <- binary_matrix[, -1]
rownames(upset_data) <- binary_matrix$Gene

# Generate UpSet plot
p1 <- upset(upset_data,
            sets = colnames(upset_data),
            nset = 4,
            order.by = "freq",
            mainbar.y.label = "No. of genes",
            sets.x.label = "No. of genes",
            decreasing = TRUE)

# Save plots
png(file.path(output_dir, "upset.png"), width = 6, height = 4, units = "in", res = 300)
print(p1)
dev.off()

pdf(file.path(output_dir, "upset.pdf"), width = 6, height = 4)
print(p1)
dev.off()

# Format and save Excel
excel_output <- binary_matrix %>%
  mutate(
    List1 = ifelse(List1 == 1, "Yes", ""),
    List2 = ifelse(List2 == 1, "Yes", ""), 
    List3 = ifelse(List3 == 1, "Yes", ""),
    List4 = ifelse(List4 == 1, "Yes", "")
  )

write.xlsx(excel_output, file.path(output_dir, "gene_lists_intersection.xlsx"))

# Print statistics
cat("Statistics:\n")
cat("Total unique genes:", nrow(binary_matrix), "\n")
cat("List1 only:", sum(binary_matrix$List1 == 1 & rowSums(binary_matrix[,2:5]) == 1), "\n")
cat("List2 only:", sum(binary_matrix$List2 == 1 & rowSums(binary_matrix[,2:5]) == 1), "\n")
cat("List3 only:", sum(binary_matrix$List3 == 1 & rowSums(binary_matrix[,2:5]) == 1), "\n")
cat("List4 only:", sum(binary_matrix$List4 == 1 & rowSums(binary_matrix[,2:5]) == 1), "\n")
cat("In all four lists:", sum(rowSums(binary_matrix[,2:5]) == 4), "\n")
cat("In three or more lists:", sum(rowSums(binary_matrix[,2:5]) >= 3), "\n")