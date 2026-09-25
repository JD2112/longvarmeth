#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(VennDiagram)
  library(ggplot2)
  library(rlang)
  library(grid)
})

compare_variants <- function(vcf1_path, vcf2_path, label1, label2, output_prefix) {
  message("\n------------------------------------------------------------")
  message(paste0("🔍 Comparing ", label1, " vs ", label2, " ..."))
  message("------------------------------------------------------------")

  #========================
  # Read both VCF files
  #========================
  vcf1 <- readVcf(vcf1_path)
  vcf2 <- readVcf(vcf2_path)

  #========================
  # Extract variant info
  #========================
  extract_vcf_df <- function(vcf) {
    data.frame(
      chr  = as.character(seqnames(rowRanges(vcf))),
      pos  = start(rowRanges(vcf)),
      ref  = as.character(ref(vcf)),
      alt  = sapply(alt(vcf), function(x) paste(x, collapse = ",")),
      qual = qual(vcf)
    )
  }

  df1 <- extract_vcf_df(vcf1)
  df2 <- extract_vcf_df(vcf2)

  #========================
  # Merge variant sets
  #========================
  merged <- full_join(
    df1, df2,
    by = c("chr", "pos", "ref", "alt"),
    suffix = c(paste0("_", label1), paste0("_", label2))
  )

  col1 <- paste0("qual_", label1)
  col2 <- paste0("qual_", label2)

  #========================
  # Determine shared/unique
  #========================
  merged$status <- case_when(
    !is.na(merged[[col1]]) & !is.na(merged[[col2]]) ~ "Shared",
    !is.na(merged[[col1]]) &  is.na(merged[[col2]]) ~ paste0("Unique_", label1),
    is.na(merged[[col1]])  & !is.na(merged[[col2]]) ~ paste0("Unique_", label2),
    TRUE ~ NA_character_
  )

  #========================
  # Summary
  #========================
  summary <- merged %>%
    count(status, name = "Variant_Count") %>%
    arrange(desc(Variant_Count))

  print(summary)

  #========================
  # Write results
  #========================
  write.csv(merged, paste0(output_prefix, "_variant_comparison.csv"), row.names = FALSE)

  #========================
  # Venn Diagram
  #========================
  area1 <- nrow(df1)
  area2 <- nrow(df2)
  cross_area <- sum(merged$status == "Shared", na.rm = TRUE)

  venn <- draw.pairwise.venn(
    area1 = area1,
    area2 = area2,
    cross.area = min(cross_area, min(area1, area2)),
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

  message(paste0("✅ Comparison done for ", label1, " vs ", label2, "."))
}


#=====================
# Run comparisons
#=====================
compare_variants(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21.vcf.gz",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22.vcf.gz",
  "21", "22", "nuclear_compare_21_22"
)

compare_variants(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24.vcf.gz",
  "23", "24", "nuclear_compare_23_24"
)

message("🎯 All nuclear variant comparisons completed successfully.")
