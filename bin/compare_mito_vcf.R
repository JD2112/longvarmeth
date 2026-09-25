#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(VennDiagram)
  library(ggplot2)
  library(rlang)
})

compare_variants <- function(vcf1_path, vcf2_path, label1, label2, output_prefix) {
  message(paste0("Comparing ", label1, " vs ", label2, " ..."))

  # Read both VCF files
  vcf1 <- readVcf(vcf1_path)
  vcf2 <- readVcf(vcf2_path)

  # Extract variant data
  df1 <- data.frame(
    chr = as.character(seqnames(rowRanges(vcf1))),
    pos = start(rowRanges(vcf1)),
    ref = as.character(ref(vcf1)),
    alt = as.character(unlist(alt(vcf1))),
    qual = qual(vcf1)
  )

  df2 <- data.frame(
    chr = as.character(seqnames(rowRanges(vcf2))),
    pos = start(rowRanges(vcf2)),
    ref = as.character(ref(vcf2)),
    alt = as.character(unlist(alt(vcf2))),
    qual = qual(vcf2)
  )

  # Merge
  merged <- full_join(
    df1, df2,
    by = c("chr", "pos", "ref", "alt"),
    suffix = c(paste0("_", label1), paste0("_", label2))
  )

  # Build dynamic column names
  col1 <- paste0("qual_", label1)
  col2 <- paste0("qual_", label2)

  # Determine shared and unique variants safely
  merged$status <- case_when(
    !is.na(merged[[col1]]) & !is.na(merged[[col2]]) ~ "Shared",
    !is.na(merged[[col1]]) & is.na(merged[[col2]]) ~ paste0("Unique_", label1),
    is.na(merged[[col1]]) & !is.na(merged[[col2]]) ~ paste0("Unique_", label2),
    TRUE ~ NA_character_
  )

  # Summary counts
  summary <- merged %>%
    count(status, name = "Variant_Count")

  print(summary)

  # Save results
  write.csv(merged, paste0(output_prefix, "_variant_comparison.csv"), row.names = FALSE)

  # Safe cross area
  cross_area <- sum(merged$status == "Shared", na.rm = TRUE)

  # Draw Venn diagram
  venn <- draw.pairwise.venn(
    area1 = nrow(df1),
    area2 = nrow(df2),
    cross.area = min(cross_area, min(nrow(df1), nrow(df2))), # ensure valid
    category = c(label1, label2),
    fill = c("#66c2a5", "#fc8d62"),
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = c(-20, 20)
  )

  png(paste0(output_prefix, "_venn.png"), width = 1200, height = 1000, res = 200)
  grid.draw(venn)
  dev.off()

  message(paste0("✔ Comparison done for ", label1, " vs ", label2, "."))
}



#=====================
# Run comparisons
#=====================
compare_variants(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/mitochondrial/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21.vcf.gz",
  "variants/mitochondrial/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22.vcf.gz",
  "21", "22", "mito_compare_21_22"
)

compare_variants(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/mitochondrial/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz",
  "variants/mitochondrial/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24.vcf.gz",
  "23", "24", "mito_compare_23_24"
)

message("✅ All mitochondrial variant comparisons completed.")
