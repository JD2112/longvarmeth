#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(VariantAnnotation)
  library(dplyr)
})

#------------------------
# Function: methylation-to-mutation link
#------------------------
meth_to_mut <- function(methyl_bed, vcf_file, label, output_prefix) {
  message("\n------------------------------------------------------------")
  message(paste0("🔬 Checking methylation→mutation overlap for ", label, " ..."))
  message("------------------------------------------------------------")

  # Load methylation sites (BED)
  meth <- import(methyl_bed, format = "BED")
  meth_gr <- GRanges(seqnames = seqnames(meth),
                     ranges = ranges(meth))

  # Load variants
  vcf <- readVcf(vcf_file)
  var_gr <- GRanges(seqnames = seqnames(rowRanges(vcf)),
                    ranges = ranges(rowRanges(vcf)))

  # Find overlaps
  overlaps <- findOverlaps(meth_gr, var_gr)
  meth_hits <- meth_gr[queryHits(overlaps)]
  var_hits  <- var_gr[subjectHits(overlaps)]

  if (length(overlaps) == 0) {
    message(paste0("⚠️ No overlaps found for ", label))
    return(NULL)
  }

  # Combine data
  combined <- data.frame(
    chr = as.character(seqnames(meth_hits)),
    start = start(meth_hits),
    end = end(meth_hits),
    variant_pos = start(var_hits),
    ref = as.character(ref(vcf)[subjectHits(overlaps)]),
    alt = sapply(alt(vcf)[subjectHits(overlaps)], function(x) paste(x, collapse = ",")),
    qual = qual(vcf)[subjectHits(overlaps)]
  )

  write.csv(combined, paste0(output_prefix, "_meth_to_mut_overlap.csv"), row.names = FALSE)

  message(paste0("✅ ", nrow(combined), " methylation sites overlap mutation positions for ", label, "."))
  invisible(combined)
}


#=====================
# Run comparisons
#=====================

meth_to_mut(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21_mods.bed",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21.vcf.gz",
  "barcode21",
  "nuclear_barcode21"
)

meth_to_mut(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22_mods.bed",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22.vcf.gz",
  "barcode22",
  "nuclear_barcode22"
)

meth_to_mut(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23_mods.bed",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz",
  "barcode23",
  "nuclear_barcode23"
)

meth_to_mut(
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24_mods.bed",
  "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24.vcf.gz",
  "barcode24",
  "nuclear_barcode24"
)

message("🎯 All methylation→mutation overlap analyses completed.")
