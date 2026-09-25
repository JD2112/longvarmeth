#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(VennDiagram)
  library(grid)
})

message("============================================================")
message("🧬 Progressive methylation → mutation analysis (mitochondrial)")
message("============================================================")

# ---------- Paths (adjust if necessary) ----------
methyl_paths_mito <- list(
  "21" = "methylation/mitochondrial_AU/*barcode21_mods.bed",
  "22" = "methylation/mitochondrial_AU/*barcode22_mods.bed",
  "23" = "methylation/mitochondrial_AU/*barcode23_mods.bed",
  "24" = "methylation/mitochondrial_AU/*barcode24_mods.bed"
)

variant_paths_mito <- list(
  "21" = "variants/mitochondrial/*barcode21.vcf.gz",
  "22" = "variants/mitochondrial/*barcode22.vcf.gz",
  "23" = "variants/mitochondrial/*barcode23.vcf.gz",
  "24" = "variants/mitochondrial/*barcode24.vcf.gz"
)

# ---------- Helpers ----------
# Fast, robust BED loader (keeps first 3 columns: chr,start,end)
load_bed_fast <- function(path, label) {
  if (!file.exists(path)) stop("Missing BED: ", path)
  df <- read_tsv(path, col_names = FALSE, comment = "#",
                 col_types = cols(.default = col_skip(),
                                  X1 = col_character(), X2 = col_double(), X3 = col_double()))
  df <- df %>% select(chr = X1, start = X2, end = X3)
  df$label <- label
  # collapse duplicates by mean start/end (just unique positions)
  df %>% distinct(chr, start, end, .keep_all = TRUE)
}

# Read VCF minimally and expand multi-ALT alleles
load_vcf_minimal <- function(path, label, genome_build = "MN240408.1") {
  if (!file.exists(path)) stop("Missing VCF: ", path)
  param <- ScanVcfParam(info = NA, geno = NA, what = c("CHROM", "POS", "REF", "ALT", "QUAL"))
  vcf <- readVcf(path, genome = genome_build, param = param)
  vr <- rowRanges(vcf)
  ref_vec <- as.character(ref(vcf))
  alt_list <- alt(vcf)

  # expand to one row per ALT allele
  rows <- lapply(seq_along(vr), function(i) {
    if (length(alt_list[[i]]) == 0) return(NULL)
    data.frame(
      chr = as.character(seqnames(vr)[i]),
      pos = start(vr)[i],
      ref = ref_vec[i],
      alt = as.character(alt_list[[i]]),
      qual = qual(vcf)[i],
      label = label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---------- Load mitochondrial methylation (all gens) ----------
message("Loading mitochondrial methylation BEDs...")
methyl_all_m <- bind_rows(lapply(names(methyl_paths_mito), function(lbl) {
  files <- Sys.glob(methyl_paths_mito[[lbl]])
  if (length(files) == 0) {
    message("  WARNING: no methylation BEDs for generation ", lbl)
    return(tibble())
  }
  bind_rows(lapply(files, load_bed_fast, label = lbl))
}))

if (nrow(methyl_all_m) == 0) {
  stop("No mitochondrial methylation data found. Check paths.")
}
message("Loaded methylation rows: ", nrow(methyl_all_m))

# Convert to GRanges (use start position as site)
methyl_gr_m <- GRanges(seqnames = methyl_all_m$chr,
                       ranges = IRanges(start = methyl_all_m$start, end = methyl_all_m$end),
                       gen = methyl_all_m$label)

# ---------- Load mitochondrial variants ----------
message("Loading mitochondrial VCFs (minimal fields)...")
mut_all_m <- bind_rows(lapply(names(variant_paths_mito), function(lbl) {
  files <- Sys.glob(variant_paths_mito[[lbl]])
  if (length(files) == 0) {
    message("  WARNING: no VCFs for generation ", lbl)
    return(data.frame())
  }
  bind_rows(lapply(files, load_vcf_minimal, label = lbl, genome_build = "MN240408.1"))
}))

if (nrow(mut_all_m) == 0) {
  message("No mitochondrial VCFs found — skipping variant overlap analysis.")
  quit(status = 0)
}
message("Loaded variant rows: ", nrow(mut_all_m))

# Convert to GRanges (variant pos as 1bp range)
mut_gr_m <- GRanges(seqnames = mut_all_m$chr,
                    ranges = IRanges(start = mut_all_m$pos, end = mut_all_m$pos),
                    ref = mut_all_m$ref,
                    alt = mut_all_m$alt,
                    qual = mut_all_m$qual,
                    gen = mut_all_m$label)

# ---------- Find overlaps: methylation site overlaps variant position ----------
message("Finding overlaps (methylation -> mutation) ...")
tol <- 0  # exact overlap at same coordinate; increase if you want +/- tolerance
ov_m <- findOverlaps(methyl_gr_m, mut_gr_m, maxgap = tol)

if (length(ov_m) == 0) {
  message("No overlaps found between methylation sites and variants in mitochondrial data.")
} else {
  df_m <- data.frame(
    chr = as.character(seqnames(methyl_gr_m[queryHits(ov_m)])),
    meth_start = start(methyl_gr_m[queryHits(ov_m)]),
    meth_end = end(methyl_gr_m[queryHits(ov_m)]),
    methyl_gen = mcols(methyl_gr_m[queryHits(ov_m)])$gen,
    var_pos = start(mut_gr_m[subjectHits(ov_m)]),
    var_ref = mcols(mut_gr_m[subjectHits(ov_m)])$ref,
    var_alt = mcols(mut_gr_m[subjectHits(ov_m)])$alt,
    var_gen = mcols(mut_gr_m[subjectHits(ov_m)])$gen,
    stringsAsFactors = FALSE
  )

  # Summarize progression per coordinate
  trend_m <- df_m %>%
    group_by(chr, var_pos, var_ref, var_alt) %>%
    summarise(
      methyl_in = paste(sort(unique(methyl_gen)), collapse = ","),
      mutate_in = paste(sort(unique(var_gen)), collapse = ","),
      n_meth = n_distinct(methyl_gen),
      n_mut = n_distinct(var_gen),
      .groups = "drop"
    )

  write.csv(df_m, "mito_meth_to_mut_overlaps_raw.csv", row.names = FALSE)
  write.csv(trend_m, "mito_progressive_meth_to_mut_summary.csv", row.names = FALSE)

  message("✅ Saved overlap table: mito_meth_to_mut_overlaps_raw.csv")
  message("✅ Saved progression summary: mito_progressive_meth_to_mut_summary.csv")

  # quick summary counts
  summary_counts <- trend_m %>%
    mutate(meth_first = as.integer(sub(",.*", "", methyl_in)),
           mut_first = as.integer(sub(",.*", "", mutate_in)),
           relation = case_when(
             mut_first > meth_first ~ "Methylation_precedes_mutation",
             mut_first == meth_first ~ "Simultaneous",
             mut_first < meth_first ~ "Mutation_precedes_methylation",
             TRUE ~ "Other"
           )) %>%
    count(relation)

  print(summary_counts)

  # Venn-like counts for a specific pair (21 vs 22) example
  # (you can change labels for other pairings)
  get_pair_venn <- function(a, b) {
    meth_a <- methyl_all_m %>% filter(label == a) %>% mutate(site = paste(chr, start, end, sep=":"))
    meth_b <- methyl_all_m %>% filter(label == b) %>% mutate(site = paste(chr, start, end, sep=":"))
    var_a <- mut_all_m %>% filter(label == a) %>% mutate(site = paste(chr, pos, pos, sep=":"))
    var_b <- mut_all_m %>% filter(label == b) %>% mutate(site = paste(chr, pos, pos, sep=":"))

    # methylation overlap between gen a and b
    m_shared <- length(intersect(meth_a$site, meth_b$site))
    m_a <- nrow(meth_a); m_b <- nrow(meth_b)

    # variant overlap between gen a and b
    v_shared <- length(intersect(var_a$site, var_b$site))
    v_a <- nrow(var_a); v_b <- nrow(var_b)

    list(meth = c(a=m_a, b=m_b, shared=m_shared),
         var = c(a=v_a, b=v_b, shared=v_shared))
  }

  pair_21_22 <- get_pair_venn("21","22")
  message("Pair 21 vs 22 methylation (a,b,shared): ", paste(pair_21_22$meth, collapse = ", "))
  message("Pair 21 vs 22 variants  (a,b,shared): ", paste(pair_21_22$var, collapse = ", "))

  # optional Venn for methylation 21 vs 22
  venn_meth <- draw.pairwise.venn(pair_21_22$meth["a"], pair_21_22$meth["b"], cross.area = pair_21_22$meth["shared"],
                                  category = c("21","22"), fill = c("lightblue","lightgreen"), alpha = c(0.5,0.5))
  png("mito_methyl_21_22_venn.png", width = 800, height = 800)
  grid.draw(venn_meth); dev.off()
  message("Saved mito_methyl_21_22_venn.png")
}

message("\n✅ Mitochondrial progressive meth→mut analysis completed.")
