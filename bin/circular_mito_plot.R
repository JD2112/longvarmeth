#==========================
# Circular Mitochondrial Genome Plot (PNG output)
#==========================

library(circlize)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

#==========================
# USER CONFIGURATION
#==========================
mito_gff_file <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/reference_genome/MN240408_1.gff3"
sites_csv <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/differential_pathway_results/mitochondrial_AU_barcode21_vs_barcode23_differential_annotated.csv"

significant_col <- "significant"
delta_col <- "delta"
gene_label_cex <- 0.6
out_png <- "mitochondrial_genome_plot.png"

#==========================
# LOAD DATA
#==========================
mito_gff <- import(mito_gff_file)
genes <- mito_gff[mito_gff$type == "gene"]
genome_length <- max(end(genes))

sites <- read.csv(sites_csv, stringsAsFactors = FALSE) %>%
  mutate(
    pos = as.numeric(pos),
    delta = as.numeric(delta),
    significant = as.logical(significant)
  ) %>%
  filter(pos >= 1 & pos <= genome_length)

#==========================
# CREATE GRANGES
#==========================
genes_gr <- GRanges(
  seqnames = rep("mt", length(genes)),
  ranges = IRanges(start = start(genes), end = end(genes)),
  gene_id = genes$gene_id,
  strand = strand(genes)
)

sites_gr <- GRanges(
  seqnames = rep("mt", nrow(sites)),
  ranges = IRanges(start = sites$pos, end = sites$pos)
)

# Distance to nearest gene (optional)
nearest_idx <- nearest(sites_gr, genes_gr)
sites$distance_to_gene <- distance(sites_gr, genes_gr[nearest_idx])

#==========================
# PLOT
#==========================
# PNG output
png(out_png, width = 2000, height = 2000, res = 300)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.01, 0.01))

# Initialize layout
circos.initialize(factors = "mt", xlim = c(0, genome_length))

# Track 1: genes
circos.trackPlotRegion(factors = "mt", ylim = c(0,1), track.height = 0.15, bg.border = NA,
                       panel.fun = function(x, y) {
  for(i in seq_along(genes)) {
    g <- genes[i]
    # Gene rectangle
    circos.rect(start(g), 0, end(g), 1,
                sector.index = "mt",
                col = ifelse(strand(g) == "+", "skyblue", "salmon"),
                border = "black")
    # Gene label
    circos.text((start(g)+end(g))/2, 0.5,
                labels = genes$Name[i],
                sector.index = "mt",
                facing = "clockwise", niceFacing = TRUE, cex = gene_label_cex)
  }
})

# Track 2: variants
circos.trackPlotRegion(factors = "mt", ylim = c(0,1), track.height = 0.25, bg.border = NA,
                       panel.fun = function(x, y) {
  max_delta <- max(abs(sites$delta), na.rm = TRUE)
  for(j in seq_len(nrow(sites))) {
    s <- sites[j, ]
    if(!is.na(s$delta) & !is.na(s$significant)) {
      col_pt <- ifelse(s$significant, "red", "gray50")
      # Point size scaled by delta magnitude
      cex_pt <- 0.5 + 0.5 * min(abs(s$delta)/max_delta, 1)
      circos.points(s$pos, 0.5, sector.index = "mt",
                    col = col_pt, pch = 16, cex = cex_pt)
    }
  }
})

# Axis
circos.axis(h = "top", major.at = seq(0, genome_length, by = 500),
            labels.cex = 0.6, labels.facing = "clockwise")

dev.off()
cat("✅ PNG plot saved to:", out_png, "\n")
