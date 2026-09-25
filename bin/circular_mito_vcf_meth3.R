library(circlize)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(VariantAnnotation)

# User config
mito_gff_file <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/reference_genome/MN240408_1.gff3"
sites_csv       <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/differential_pathway_results/mitochondrial_AU_barcode21_vs_barcode23_differential_annotated.csv"
vcf_file        <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/mitochondrial/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz"

out_png <- "mitochondrial_genome_plot_circlize_style2.png"
gene_label_cex <- 0.6

# Load data
mito_gff <- import(mito_gff_file)
genes <- mito_gff[mito_gff$type == "gene"]
genome_length <- max(end(genes))

sites <- read.csv(sites_csv, stringsAsFactors = FALSE) %>%
  mutate(pos          = as.numeric(pos),
         delta        = as.numeric(delta),
         significant  = as.logical(significant)) %>%
  filter(pos >= 1 & pos <= genome_length)

vcf <- readVcf(vcf_file)
vcf_gr <- rowRanges(vcf)
ref_alleles <- as.character(ref(vcf))
alt_alleles <- as.character(unlist(alt(vcf)))
vcf_df <- data.frame(
  pos = start(vcf_gr),
  ref = ref_alleles,
  alt = alt_alleles,
  stringsAsFactors = FALSE
) %>%
  mutate(is_snp = (nchar(ref) == 1 & nchar(alt) == 1),
         variant_type = ifelse(is_snp, "SNP", "Indel"))

# Prepare dataframes
genes_df <- data.frame(
  chr       = "mt",
  start     = start(genes),
  end       = end(genes),
  strand    = as.character(strand(genes)),
  gene_name = genes$Name,
  stringsAsFactors = FALSE
)

delta_df <- sites %>%
  transmute(chr = "mt",
            start = pos,
            end   = pos,
            delta = delta,
            significant = significant)

vcf_snp_df   <- vcf_df %>% filter(variant_type == "SNP")   %>% transmute(chr="mt", start=pos, end=pos)
vcf_indel_df<- vcf_df %>% filter(variant_type == "Indel") %>% transmute(chr="mt", start=pos, end=pos)

# Plot
# ==========================
# Circular Mitochondrial Genome Plot with Gene Bars & Variants
# ==========================

png(out_png, width = 2500, height = 2500, res = 300)
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree   = 1,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.01, 0.01)  # smaller top/bottom margins for each track
)

# Initialize
circos.initialize(factors = "mt", xlim = c(0, genome_length))

# Track 1: Axis
circos.trackPlotRegion(factors="mt",
                       ylim = c(0,1),
                       track.height = 0.04,
                       bg.border   = "black",
                       panel.fun = function(x, y) {
                         major.at   <- seq(0, genome_length, by = 2000)
                         circos.axis(h = "top", major.at = major.at,
                                     labels = major.at, labels.cex = 0.5,
                                     labels.facing = "clockwise", col="gray40")
                       })

# Track 2: Genes as rectangles
circos.trackPlotRegion(factors="mt",
                       ylim   = c(0,1),
                       track.height = 0.10,
                       bg.border   = NA,
                       panel.fun = function(x, y) {
                         for(i in seq_len(nrow(genes_df))) {
                           g <- genes_df[i, ]
                           col_fill <- ifelse(g$strand == "+", "skyblue", "salmon")
                           circos.rect(g$start, 0, g$end, 1,
                                       sector.index = "mt",
                                       col = col_fill, border="black")
                         }
                       })

# Track 3: Delta sites as bars
# Track 3: Delta sites as bars (fixed)
max_delta <- max(abs(delta_df$delta), na.rm = TRUE)

circos.trackPlotRegion(
  factors = "mt",
  ylim   = c(-max_delta, max_delta),
  track.height = 0.18,
  bg.border   = NA,
  panel.fun = function(x, y) {
    df <- delta_df
    df <- df[order(df$delta), ]
    col_vec <- ifelse(df$significant, "red", "gray60")
    
    # draw bars manually with rects (more flexible than circos.barplot)
    for (j in seq_len(nrow(df))) {
      s <- df[j, ]
      circos.rect(xleft = s$start - 5,  # half-width for better visibility
                  ybottom = 0,
                  xright = s$start + 5,
                  ytop = s$delta,
                  sector.index = "mt",
                  col = col_vec[j],
                  border = NA)
    }
    
    # horizontal baseline
    circos.lines(c(0, genome_length), c(0, 0), col = "black", lwd = 0.6)
  }
)

# Track 4: Variants
# Track 4: Variants (fixed ylim)
circos.trackPlotRegion(
  factors = "mt",
  ylim        = c(0, 1.2),
  track.height = 0.06,
  bg.border    = NA,
  panel.fun = function(x, y) {
    if (nrow(vcf_snp_df) > 0)
      circos.points(vcf_snp_df$start, 0.8, col = "darkgreen", pch = 17, cex = 0.6)
    if (nrow(vcf_indel_df) > 0)
      circos.points(vcf_indel_df$start, 0.4, col = "darkorchid", pch = 15, cex = 0.6)
  }
)

# Track 5: Gene labels outside (improved)
# Track 5: Gene labels outside (larger, visible)
circos.trackPlotRegion(
  factors = "mt",
  ylim         = c(0, 1.5),   # more space for labels
  track.height = 0.10,        # slightly taller
  bg.border    = NA,
  panel.fun    = function(x, y) {
    for (i in seq_len(nrow(genes_df))) {
      g <- genes_df[i, ]
      mid <- (g$start + g$end) / 2
      offset <- 1.1 + 0.2 * (i %% 2)  # alternate outward offset
      circos.text(
        mid, offset,
        labels = g$gene_name,
        sector.index = "mt",
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.9,        # 🔹 larger font size
        col = "black",
        font = 2          # bold
      )
    }
  }
)

# Legend
legend("topright",
       legend = c("Gene (+ strand)", "Gene (− strand)",
                  "Significant Δ", "Non-significant Δ", "SNP", "Indel"),
       pch    = c(22,22,15,15,17,15),
       pt.bg  = c("skyblue","salmon","red","gray60","darkgreen","darkorchid"),
       col    = c("black","black","red","gray60","darkgreen","darkorchid"),
       cex    = 0.8, bty="n",
       title  = "Legend")

dev.off()
cat("✅ Circos plot saved to:", out_png, "\n")
