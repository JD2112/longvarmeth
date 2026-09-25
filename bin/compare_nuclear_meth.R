#!/usr/bin/env Rscript
# ===================================================================
# Compare nuclear methylation: 21 vs 22, 23 vs 24
# Author: Jyotirmoy Das
# Date: 2025-10-27
# ===================================================================

library(tidyverse)
library(VennDiagram)
library(ggrepel)
library(grid)

# ============ Helper function ============

compare_methylation <- function(file1, file2, label1, label2, out_prefix) {
  # Read the BED files
  df1 <- read_tsv(file1, col_names = FALSE, show_col_types = FALSE)
  df2 <- read_tsv(file2, col_names = FALSE, show_col_types = FALSE)
  
  # Select relevant columns
  df1 <- df1 %>% select(chr = X1, start = X2, end = X3, strand = X6, mod_fraction = X5)
  df2 <- df2 %>% select(chr = X1, start = X2, end = X3, strand = X6, mod_fraction = X5)
  
  # Collapse duplicates (average methylation per CpG)
  df1 <- df1 %>%
    group_by(chr, start, end, strand) %>%
    summarise(mod_fraction = mean(mod_fraction, na.rm = TRUE), .groups = "drop")
  
  df2 <- df2 %>%
    group_by(chr, start, end, strand) %>%
    summarise(mod_fraction = mean(mod_fraction, na.rm = TRUE), .groups = "drop")
  
  # Merge unique CpGs
  merged <- inner_join(df1, df2, by = c("chr", "start", "end", "strand"),
                       suffix = c(paste0("_", label1), paste0("_", label2)))
  
  # Compute delta
  merged <- merged %>%
    mutate(delta = .data[[paste0("mod_fraction_", label2)]] -
                   .data[[paste0("mod_fraction_", label1)]])
  
  # Save merged table
  write_tsv(merged, paste0(out_prefix, "_merged.tsv"))
  
  # ---------- 1. Histogram ----------
  p1 <- ggplot(merged, aes(x = delta)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Delta Methylation (", label2, " vs ", label1, ")"),
         x = paste0("Change in methylation fraction (", label2, " - ", label1, ")"),
         y = "CpG count")
  
  ggsave(paste0(out_prefix, "_delta_hist.png"), p1, width = 7, height = 5, dpi = 300)
  
  # ---------- 2. Scatter plot ----------
  p2 <- ggplot(merged, aes(
    x = .data[[paste0("mod_fraction_", label1)]],
    y = .data[[paste0("mod_fraction_", label2)]]
  )) +
    geom_point(alpha = 0.4) +
    geom_abline(linetype = "dashed", color = "red") +
    theme_classic(base_size = 14) +
    labs(title = paste0("CpG Methylation Concordance (", label1, " vs ", label2, ")"),
         x = paste0("Generation ", label1),
         y = paste0("Generation ", label2))
  
  ggsave(paste0(out_prefix, "_scatter.png"), p2, width = 7, height = 5, dpi = 300)
  
  # ---------- 3. Venn Diagram ----------
  venn.plot <- draw.pairwise.venn(
    area1 = nrow(df1),
    area2 = nrow(df2),
    cross.area = nrow(merged),
    category = c(paste("Barcode", label1), paste("Barcode", label2)),
    fill = c("lightblue", "lightgreen"),
    alpha = c(0.5, 0.5),
    cat.pos = c(0, 0)
  )
  
  png(paste0(out_prefix, "_venn.png"), width = 800, height = 800)
  grid.draw(venn.plot)
  dev.off()
  
  message("✅ Comparison completed: ", label1, " vs ", label2)
}



# ===================================================================
# Input files
# ===================================================================
meth_files <- list(
  "21" = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21_mods.bed",
  "22" = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22_mods.bed",
  "23" = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23_mods.bed",
  "24" = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24_mods.bed"
)

# ===================================================================
# Run comparisons
# ===================================================================

compare_methylation(meth_files$`21`, meth_files$`22`, "21", "22", "compare_nuclear_21vs22")
compare_methylation(meth_files$`23`, meth_files$`24`, "23", "24", "compare_nuclear_23vs24")

# ===================================================================
# 4. Summary plot: Average methylation per barcode
# ===================================================================

all_meth <- tibble(
  sample = c("21", "22", "23", "24"),
  mean_methyl = sapply(meth_files, function(f) {
    df <- read_tsv(f, col_names = FALSE, show_col_types = FALSE)
    mean(df$X5, na.rm = TRUE)
  })
)

p_summary <- ggplot(all_meth, aes(x = sample, y = mean_methyl, group = 1)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "darkred") +
  theme_minimal(base_size = 14) +
  labs(title = "Average nuclear methylation across generations",
       x = "Generation (Barcode)", y = "Mean methylation fraction")

ggsave("compare_nuclear_mean_methylation.png", p_summary, width = 7, height = 5, dpi = 300)

message("✅ All comparisons and summary plot completed.")
