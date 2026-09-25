library(circlize)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(VariantAnnotation)

# User config
mito_gff_file <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/reference_genome/MN240408_1.gff3"
sites_csv       <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/differential_pathway_results/mitochondrial_AU_barcode21_vs_barcode23_differential_annotated.csv"
vcf_file        <- "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/variants/mitochondrial/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23.vcf.gz"

out_png <- "mitochondrial_genome_plot_circlize_style.png"
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
png(out_png, width = 2500, height = 2500, res = 300)
circos.clear()
circos.par(start.degree = 90,
           gap.degree   = 1,
           cell.padding = c(0,0,0,0))

# Initialize sectors
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

# Track 3: Delta sites (many points)
max_delta <- max(abs(delta_df$delta), na.rm = TRUE)
circos.trackPlotRegion(factors = "mt",
                       ylim   = c(0,1),
                       track.height = 0.12,
                       bg.border   = NA,
                       panel.fun = function(x, y) {
                         for(j in seq_len(nrow(delta_df))) {
                           s <- delta_df[j, ]
                           col_pt <- ifelse(s$significant, "red", "gray50")
                           cex_pt <- 0.4 + 0.6 * min(abs(s$delta)/max_delta, 1)
                           circos.points(s$start, 0.5, col = col_pt, pch = 16, cex = cex_pt)
                         }
                       })

# Track 4: VCF Variants (SNP + Indel)
circos.trackPlotRegion(factors="mt",
                       ylim        = c(0,1),
                       track.height = 0.06,
                       bg.border    = NA,
                       panel.fun = function(x, y) {
                         if(nrow(vcf_snp_df)   > 0) circos.points(vcf_snp_df$start,   0.6, col="darkgreen", pch=17, cex=0.5)
                         if(nrow(vcf_indel_df)> 0) circos.points(vcf_indel_df$start,0.4, col="darkorchid", pch=15, cex=0.5)
                       })

# Track 5: Gene labels outside
circos.trackPlotRegion(factors="mt",
                       ylim         = c(0,1),
                       track.height = 0.08,
                       bg.border    = NA,
                       panel.fun    = function(x, y) {
                         for(i in seq_len(nrow(genes_df))) {
                           g <- genes_df[i, ]
                           mid <- (g$start + g$end)/2
                           circos.text(mid, 1.1, labels = g$gene_name,
                                       sector.index = "mt",
                                       facing = "clockwise",
                                       niceFacing = TRUE,
                                       cex = gene_label_cex)
                         }
                       })

# Legend
legend("topright",
       legend = c("Gene (+ strand)", "Gene (− strand)", "Significant delta", "Non-significant delta", "SNP", "Indel"),
       pch    = c(22,22,16,16,17,15),
       pt.bg  = c("skyblue","salmon","red","gray50","darkgreen","darkorchid"),
       col    = c("black","black","red","gray50","darkgreen","darkorchid"),
       cex    = 0.8, bty="n",
       title  = "Legend")

dev.off()
cat("✅ Plot saved to:", out_png, "\n")