library(circlize)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(VariantAnnotation)

# User config
mito_gff_file <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/reference_genome/MN240408_1.gff3"
sites_csv       <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/differential_pathway_results/mitochondrial_AU_barcode21_vs_barcode23_differential_annotated.csv"
vcf_file        <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/mitochondrial/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz"

out_png <- "mito_genome_plot_styleTWGBS.png"
gene_label_cex <- 0.6

# --- READ GFF ---
mito_gff <- import(mito_gff_file)
genes <- mito_gff[mito_gff$type == "gene"]
genome_length <- max(end(genes))

# --- READ DELTA SITES ---
sites <- read.csv(sites_csv, stringsAsFactors = FALSE) %>%
  mutate(pos = as.numeric(pos),
         delta = as.numeric(delta),
         significant = as.logical(significant)) %>%
  filter(!is.na(pos) & pos >= 1 & pos <= genome_length)

# --- READ VCF ---
vcf <- readVcf(vcf_file)
vcf_gr <- rowRanges(vcf)
ref_alleles <- as.character(ref(vcf))
alt_alleles <- as.character(unlist(alt(vcf)))
vcf_df <- data.frame(pos = start(vcf_gr),
                     ref = ref_alleles,
                     alt = alt_alleles,
                     stringsAsFactors = FALSE) %>%
  mutate(is_snp = (nchar(ref) == 1 & nchar(alt) == 1),
         variant_type = ifelse(is_snp, "SNP", "Indel"))

# --- PREPARE DATAFRAMES ---
genes_df <- data.frame(chr = "mt",
                       start = start(genes),
                       end   = end(genes),
                       strand = as.character(strand(genes)),
                       gene_name = genes$Name)

delta_df <- sites %>%
  transmute(chr = "mt", start = pos, end = pos,
            delta = delta, significant = significant)

vcf_snp_df   <- vcf_df %>% filter(variant_type == "SNP")   %>%
  transmute(chr = "mt", start = pos, end = pos)
vcf_indel_df <- vcf_df %>% filter(variant_type == "Indel") %>%
  transmute(chr = "mt", start = pos, end = pos)

#==========================
# PLOT
#==========================
png(out_png, width = 2400, height = 2400, res = 300)
circos.clear()
circos.par(start.degree = 90, gap.degree = 2,
           cell.padding = c(0, 0, 0, 0), track.margin = c(0.002, 0.002))

# Initialize for mt genome
circos.initialize(factors = "mt", xlim = c(0, genome_length))

# Track 1: base axis
circos.trackPlotRegion(factors = "mt", ylim = c(0, 1),
  track.height = 0.04, bg.border = "black",
  panel.fun = function(x, y) {
    major.at <- seq(0, genome_length, by = round(genome_length / 10))
    circos.axis(h = "top", major.at = major.at, labels = major.at,
                labels.cex = 0.5, labels.facing = "clockwise", col = "gray40")
  })

# Track 2: genes
circos.genomicTrack(genes_df, ylim = c(0, 1), track.height = 0.10, bg.border = NA,
  panel.fun = function(region, value, ...) {
    if (nrow(value) > 0) {
      for (i in seq_len(nrow(value))) {
        g <- value[i, ]
        col_fill <- ifelse(g$strand == "+", "skyblue", "salmon")
        circos.genomicRect(region[i, , drop = FALSE], value[i, , drop = FALSE],
                           col = col_fill, border = "black", lwd = 0.4)
      }
    }
  })

# Track 3: delta sites (TWGBS-DMR style)
max_delta <- max(abs(delta_df$delta), na.rm = TRUE)
circos.genomicTrack(delta_df, ylim = c(-max_delta, max_delta),
  track.height = 0.12, bg.col = "#f7f7f7", bg.border = NA,
  panel.fun = function(region, value, ...) {
    if (nrow(value) > 0) {
      h_vals <- seq(-max_delta, max_delta, by = max_delta / 4)
      for (h in h_vals)
        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#777777")
      for (i in seq_len(nrow(value))) {
        v <- value[i, ]
        col_pt <- ifelse(v$significant, "#E41A1C", "#377EB8")
        cex_pt <- 0.4 + 0.6 * min(abs(v$delta) / max_delta, 1)
        circos.points(v$start, v$delta, col = col_pt, pch = 16, cex = cex_pt)
      }
    }
  })

# Track 4: variants
circos.genomicTrack(vcf_snp_df, ylim = c(0, 1),
  track.height = 0.06, bg.border = NA,
  panel.fun = function(region, value, ...) {
    if (nrow(vcf_snp_df) > 0)
      circos.points(vcf_snp_df$start, rep(0.6, nrow(vcf_snp_df)),
                    col = "darkgreen", pch = 17, cex = 0.5)
    if (nrow(vcf_indel_df) > 0)
      circos.points(vcf_indel_df$start, rep(0.4, nrow(vcf_indel_df)),
                    col = "darkorchid", pch = 15, cex = 0.5)
  })

# Track 5: gene labels outside
circos.trackPlotRegion(factors = "mt", ylim = c(0, 1), track.height = 0.10, bg.border = NA,
  panel.fun = function(x, y) {
    for (i in seq_len(nrow(genes_df))) {
      g <- genes_df[i, ]
      mid <- (g$start + g$end) / 2
      circos.text(mid, 1.1, g$gene_name, sector.index = "mt",
                  facing = "clockwise", niceFacing = TRUE, cex = 0.6)
    }
  })

# Legend
legend("topright",
       legend = c("Gene (+)", "Gene (-)", "Δ-methylation significant", "Δ-methylation non-sig",
                  "SNP", "Indel"),
       pch = c(22, 22, 16, 16, 17, 15),
       pt.bg = c("skyblue", "salmon", "#E41A1C", "#377EB8", "darkgreen", "darkorchid"),
       col = c("black", "black", "#E41A1C", "#377EB8", "darkgreen", "darkorchid"),
       cex = 0.8, bty = "n", title = "Legend")

dev.off()
cat("✅ Plot saved to:", out_png, "\n")