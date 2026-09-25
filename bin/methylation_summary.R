#!/usr/bin/env Rscript

# Rscript scripts/methylation_summary.R
library(data.table); library(ggplot2); library(dplyr)
res_dir <- "results"
meth_dirs <- list.files("results/methylation", full.names=TRUE)
summary_list <- list()
for (gdir in meth_dirs) {
  files <- list.files(gdir, pattern="*_mods.bed", full.names=TRUE)
  if (length(files)==0) next
  for (f in files) {
    dt <- fread(f, header=FALSE)
    # expect bedMethyl-like: chrom start end name score strand ... (score at col 5)
    if (ncol(dt) < 5) next
    score <- as.numeric(dt[[5]])
    summary_list[[length(summary_list)+1]] <- data.frame(sample=basename(f), genome=basename(gdir), mean_score=mean(score, na.rm=TRUE))
  }
}
if (length(summary_list)>0) {
  df <- bind_rows(summary_list)
  p <- ggplot(df, aes(x=sample, y=mean_score, fill=genome)) + geom_col(position="dodge") + coord_flip() + theme_bw()
  dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)
  ggsave(filename="results/plots/methylation_summary_all_genomes.png", plot=p, width=12, height=6)
}