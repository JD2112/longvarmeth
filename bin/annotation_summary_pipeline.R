#!/usr/bin/env Rscript

###############################################################################
# Title: Multi-sample Annotation Summary + GO Analysis Pipeline
# Author: Jyotirmoy Das
# Description:
#   - Processes multiple methylation and variant files with shared annotation
#   - Computes distances to nearest genes, variant counts, and per-gene summaries
#   - Automatically detects organism from GFF header
#   - Performs GO enrichment on combined significant genes
###############################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(VariantAnnotation)
  library(dplyr)
  library(clusterProfiler)
  library(stringr)
  library(tools)
  library(purrr)
})

# ------------------- USER CONFIG -------------------
# Directory paths
methyl_dir <- "methylation/nuclear"
vcf_dir    <- "variants"
gff_file   <- "annotation/genes.gff3"
output_dir <- "annotation_summary_results"
dir.create(output_dir, showWarnings = FALSE)

# ------------------- ORGANISM DETECTION -------------------
cat("Detecting organism from GFF header...\n")
gff_head <- readLines(gff_file, n = 50)
organism <- "unknown"

if (any(grepl("Homo sapiens", gff_head, ignore.case = TRUE))) {
  organism <- "human"; library(org.Hs.eg.db)
} else if (any(grepl("Mus musculus", gff_head, ignore.case = TRUE))) {
  organism <- "mouse"; library(org.Mm.eg.db)
} else if (any(grepl("Danio rerio", gff_head, ignore.case = TRUE))) {
  organism <- "zebrafish"; library(org.Dr.eg.db)
} else if (any(grepl("Rattus norvegicus", gff_head, ignore.case = TRUE))) {
  organism <- "rat"; library(org.Rn.eg.db)
} else if (any(grepl("Drosophila melanogaster", gff_head, ignore.case = TRUE))) {
  organism <- "fly"; library(org.Dm.eg.db)
} else if (any(grepl("Saccharomyces cerevisiae", gff_head, ignore.case = TRUE))) {
  organism <- "yeast"; library(org.Sc.sgd.db)
} else {
  stop("❌ Could not detect organism from GFF header — please specify manually.")
}
cat(paste0("✅ Detected organism: ", organism, "\n"))

# ------------------- LOAD ANNOTATION -------------------
cat("Loading gene annotations...\n")
gff <- import(gff_file)
gene_gr <- gff[gff$type == "gene"]
gene_gr$gene_id <- ifelse(!is.na(gene_gr$gene_id),
                          gene_gr$gene_id,
                          gff$ID[gff$type == "gene"])

# ------------------- FILE LISTING -------------------
methyl_files <- list.files(methyl_dir, pattern = "_mods\\.bed$", full.names = TRUE)
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)

cat("Found methylation files:\n")
print(basename(methyl_files))
cat("Found VCF files:\n")
print(basename(vcf_files))

# Match sample names by prefix (before first underscore)
samples <- unique(str_extract(basename(methyl_files), "^[^_]+"))

# ------------------- PROCESS EACH SAMPLE -------------------
all_gene_summaries <- list()

for (sample in samples) {
  cat(paste0("\n🔹 Processing sample: ", sample, "\n"))
  
  meth_file <- methyl_files[grepl(sample, methyl_files)]
  vcf_file  <- vcf_files[grepl(sample, vcf_files)]
  
  if (length(meth_file) == 0 | length(vcf_file) == 0) {
    cat("⚠️ Skipping ", sample, " (missing BED or VCF)\n")
    next
  }
  
  # --- Load methylation data ---
  meth_df <- read.table(meth_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(meth_df)[1:3] <- c("chr", "start", "end")
  if (ncol(meth_df) >= 4) colnames(meth_df)[4] <- "value"
  
  meth_gr <- GRanges(seqnames = meth_df$chr,
                     ranges = IRanges(meth_df$start, meth_df$end),
                     value = meth_df$value)
  
  # --- Compute nearest gene ---
  nearest_idx <- nearest(meth_gr, gene_gr)
  meth_gr$nearest_gene <- gene_gr$gene_id[nearest_idx]
  meth_gr$distance_to_gene <- distance(meth_gr, gene_gr[nearest_idx])
  
  # --- Variant counts ---
  vcf <- readVcf(vcf_file)
  vcf_gr <- rowRanges(vcf)
  
  hits <- findOverlaps(vcf_gr, gene_gr)
  variant_gene_df <- data.frame(
    gene_id = gene_gr$gene_id[subjectHits(hits)],
    variant_id = names(vcf_gr[queryHits(hits)])
  )
  variant_counts <- variant_gene_df %>%
    group_by(gene_id) %>%
    summarise(variant_count = n(), .groups = "drop")
  
  # --- Per-gene summary ---
  meth_df2 <- as.data.frame(meth_gr) %>%
    select(seqnames, start, end, value, nearest_gene, distance_to_gene) %>%
    rename(chr = seqnames)
  
  gene_summary <- meth_df2 %>%
    group_by(nearest_gene) %>%
    summarise(
      mean_methylation = mean(value, na.rm = TRUE),
      n_sites = n(),
      min_distance = min(distance_to_gene, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(variant_counts, by = c("nearest_gene" = "gene_id")) %>%
    mutate(variant_count = replace_na(variant_count, 0),
           sample = sample)
  
  out_file <- file.path(output_dir, paste0(sample, "_gene_summary.tsv"))
  write.table(gene_summary, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  ✅ Saved:", out_file, "\n")
  
  all_gene_summaries[[sample]] <- gene_summary
}

# ------------------- COMBINE + GLOBAL GO -------------------
cat("\nCombining per-sample summaries...\n")

combined_summary <- bind_rows(all_gene_summaries)
combined_file <- file.path(output_dir, "combined_gene_summary.tsv")
write.table(combined_summary, combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ✅ Saved combined summary:", combined_file, "\n")

# Select top genes by variability or high methylation deviation
top_genes <- combined_summary %>%
  group_by(nearest_gene) %>%
  summarise(
    avg_methylation = mean(mean_methylation, na.rm = TRUE),
    total_variants = sum(variant_count, na.rm = TRUE)
  ) %>%
  arrange(desc(abs(avg_methylation))) %>%
  slice_head(n = 500) %>%
  pull(nearest_gene)

cat("Performing GO enrichment on top 500 genes...\n")

# Convert to Entrez IDs
gene_map <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = eval(parse(text = paste0("org.", toupper(substr(organism,1,1)), substr(organism,2,100), ".eg.db"))))

ego <- enrichGO(
  gene         = gene_map$ENTREZID,
  OrgDb        = eval(parse(text = paste0("org.", toupper(substr(organism,1,1)), substr(organism,2,100), ".eg.db"))),
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

# Save outputs
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  go_file <- file.path(output_dir, "combined_GO_enrichment.tsv")
  pdf_file <- file.path(output_dir, "combined_GO_dotplot.pdf")
  
  write.table(as.data.frame(ego), go_file, sep = "\t", row.names = FALSE, quote = FALSE)
  pdf(pdf_file, width = 8, height = 6)
  print(dotplot(ego, showCategory = 20))
  dev.off()
  
  cat("  ✅ GO enrichment and dotplot saved.\n")
}

cat("\n🎯 All analyses completed successfully!\n")
