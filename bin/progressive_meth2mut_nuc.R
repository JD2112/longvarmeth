#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(parallel)
})

message("============================================================")
message("🧬 Progressive methylation → mutation analysis (nuclear, Artemia/custom genome)")
message("============================================================")

#------------------------------------------------------------
# Paths
#------------------------------------------------------------
meth_dir <- "methylation/nuclear"
vcf_dir  <- "variants/nuclear"
output_file <- "nuclear_meth2mut.tsv"

#------------------------------------------------------------
# Load methylation BEDs safely (only barcode21–24)
#------------------------------------------------------------
message("Loading nuclear methylation BEDs (barcodes 21–24 only)...")

meth_files <- list.files(
  meth_dir,
  pattern = "barcode(21|22|23|24)_mods\\.bed$",
  full.names = TRUE
)

if (length(meth_files) == 0) stop("❌ No methylation BED files for barcodes 21–24 found in ", meth_dir)

read_bed_safely <- function(f) {
  df <- read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE)
  if (ncol(df) < 3) {
    stop("File ", f, " does not have at least 3 BED columns.")
  }
  colnames(df)[1:3] <- c("seqnames", "start", "end")
  df <- df %>% filter(!is.na(seqnames) & !is.na(start) & !is.na(end))
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  gr$barcode <- sub(".*barcode(\\d+)_mods\\.bed$", "barcode\\1", f)
  gr
}

meth_list <- lapply(meth_files, read_bed_safely)
meth <- do.call(c, meth_list)
message("Loaded methylation rows: ", length(meth))
message("Files processed:"); print(basename(meth_files))

#------------------------------------------------------------
# Load VCFs
#------------------------------------------------------------
message("Loading nuclear VCFs (minimal fields)...")
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
if (length(vcf_files) == 0) stop("❌ No VCF files found in ", vcf_dir)

vcf_list <- mclapply(vcf_files, function(f) {
  vcf <- readVcf(f)
  as(vcf, "VRanges")
}, mc.cores = max(1, parallel::detectCores() - 1))

vars <- do.call(c, vcf_list)
message("Loaded variant rows: ", length(vars))

#------------------------------------------------------------
# Harmonize sequence names
#------------------------------------------------------------
message("Checking sequence name compatibility...")
common_seq <- intersect(seqlevels(meth), seqlevels(vars))

if (length(common_seq) == 0) {
  message("⚠️ No common sequence names found — trying to standardize...")
  seqlevels(meth) <- gsub("^chr", "", seqlevels(meth))
  seqlevels(vars) <- gsub("^chr", "", seqlevels(vars))
  common_seq <- intersect(seqlevels(meth), seqlevels(vars))
}

if (length(common_seq) == 0) {
  stop("❌ Still no common sequences between methylation and variant data. 
       Check naming (e.g., contig IDs) in BED and VCF files.")
} else {
  meth <- keepSeqlevels(meth, common_seq, pruning.mode = "coarse")
  vars <- keepSeqlevels(vars, common_seq, pruning.mode = "coarse")
  message("✅ Found ", length(common_seq), " shared scaffolds/contigs.")
}

#------------------------------------------------------------
# Overlap detection
#------------------------------------------------------------
message("Finding overlaps (methylation → mutation)...")
hits <- findOverlaps(meth, vars)
message("Overlaps found: ", length(hits))

if (length(hits) > 0) {
  meth_hits <- meth[queryHits(hits)]
  var_hits  <- vars[subjectHits(hits)]
  
  # Extract safely with defaults
  get_meta <- function(x, name) {
    if (is.null(mcols(x)[[name]])) return(rep(NA, length(x)))
    v <- mcols(x)[[name]]
    if (is.list(v)) v <- sapply(v, function(i) paste(i, collapse = ","))
    v
  }
  
  df <- data.frame(
    seqname    = as.character(seqnames(meth_hits)),
    meth_start = start(meth_hits),
    meth_end   = end(meth_hits),
    barcode    = meth_hits$barcode,
    var_pos    = start(var_hits),
    var_ref    = get_meta(var_hits, "REF"),
    var_alt    = get_meta(var_hits, "ALT"),
    var_qual   = get_meta(var_hits, "QUAL")
  )
  
  # Sanity check
  stopifnot(nrow(df) == length(hits))
  
  write_tsv(df, output_file)
  message("✅ Overlaps written to: ", output_file)
  
  # Summary by barcode
  message("\n📊 Overlap summary per barcode:")
  print(df %>% count(barcode, name = "n_overlaps"))
  
} else {
  message("No overlaps found between methylation sites and variants in nuclear data.")
}

message("✅ Nuclear progressive meth→mut analysis (Artemia, barcodes 21–24) completed.")
