#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler
import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


# -----------------------------
# CONFIG
# -----------------------------
regions = {
    "nuclear": {
        "barcode21": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21_mods.bed",
        "barcode22": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22_mods.bed",
        "barcode23": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23_mods.bed",
        "barcode24": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/nuclear/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24_mods.bed",
    },
    "mitochondrial_AU": {
        "barcode21": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/mitochondrial_AU/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode21_mods.bed",
        "barcode22": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/mitochondrial_AU/92b1c360-ee89-4a61-af36-8541e3ee5415_SQK-NBD114-24_barcode22_mods.bed",
        "barcode23": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/mitochondrial_AU/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode23_mods.bed",
        "barcode24": "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results/methylation/mitochondrial_AU/c1405676-fa9d-4b3b-8ce4-2d2a1ac64627_SQK-NBD114-24_barcode24_mods.bed",
    },
}

gff_file = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/reference_genome/GCF_032884065.1_ASM3288406v1_genomic.gff"
outdir = "differential_pathway_results"
os.makedirs(outdir, exist_ok=True)


# -----------------------------
# FUNCTIONS
# -----------------------------
def read_bed(path):
    """Read BED-like modkit file"""
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, :5]  # Keep first 5 columns
    df.columns = ["chrom", "start", "end", "context", "prob"]
    df["prob"] = pd.to_numeric(df["prob"], errors="coerce")
    df = df.dropna(subset=["prob"])
    df["pos"] = df["start"]
    return df


def summarize_methylation(df):
    """Summarize average methylation probability per CpG"""
    return df.groupby(["chrom", "pos"]).agg(
        mean_prob=("prob", "mean"), n=("prob", "count")
    ).reset_index()


def differential_analysis(df1, df2, label1, label2):
    """Compute Δmethylation and p-value proxy (abs diff > threshold)"""
    merged = pd.merge(df1, df2, on=["chrom", "pos"], suffixes=(f"_{label1}", f"_{label2}"))
    merged["delta"] = merged[f"mean_prob_{label2}"] - merged[f"mean_prob_{label1}"]
    merged["significant"] = np.abs(merged["delta"]) > 0.2
    return merged


def plot_volcano(df, label1, label2, prefix):
    plt.figure(figsize=(6, 5))
    plt.scatter(df["delta"], -np.log10(df["n_" + label1] + 1), c=df["significant"].map({True: "red", False: "gray"}), alpha=0.7)
    plt.xlabel("Δ Methylation (β diff)")
    plt.ylabel("-log10(count + 1)")
    plt.title(f"Differential Methylation: {label2} vs {label1}")
    plt.tight_layout()
    plt.savefig(f"{prefix}_volcano.png", dpi=300)
    plt.close()


# -----------------------------
# MAIN PIPELINE
# -----------------------------
gp = GProfiler(return_dataframe=True)

for region, samples in regions.items():
    print(f"\n🧬 Running region: {region}")

    # Read and summarize
    dfs = {s: summarize_methylation(read_bed(p)) for s, p in samples.items()}

    comparisons = [("barcode21", "barcode23"), ("barcode22", "barcode24")]
    summary_list = []

    for a, b in comparisons:
        diff = differential_analysis(dfs[a], dfs[b], a, b)
        diff_sig = diff[diff["significant"]]

        prefix = os.path.join(outdir, f"{region}_{a}_vs_{b}")
        diff.to_csv(f"{prefix}_differential.csv", index=False)

        if not diff_sig.empty:
            plot_volcano(diff, a, b, prefix)

            # Map to genes (optional)
            annots = []
            with open(gff_file) as fh:
                for line in fh:
                    if "\tgene\t" in line:
                        parts = line.split("\t")
                        chrom, start, end = parts[0], int(parts[3]), int(parts[4])
                        gene_id = parts[-1].split(";")[0].replace("ID=", "").strip()
                        annots.append((chrom, start, end, gene_id))
            ann_df = pd.DataFrame(annots, columns=["chrom", "start", "end", "gene_id"])

            # Overlap by position
            genes_hit = pd.merge(
                diff_sig, ann_df, on="chrom", how="inner"
            )
            genes_hit = genes_hit[
                (genes_hit["pos"] >= genes_hit["start"]) & (genes_hit["pos"] <= genes_hit["end"])
            ]

            # --- After annotation step ---
            print("\n🧩 Annotating differentially methylated regions...")
            # ann_df = annotate_regions(
            #     dmr_df,
            #     annotation_file=args.annotation,
            #     genome=args.genome,
            # )

            # --- Debug: show annotation summary ---
            print(f"✅ Annotation dataframe shape: {ann_df.shape}")
            print("🧬 Annotation columns:", list(ann_df.columns))
            print(ann_df.head(10))

            # --- Extract genes overlapping the significant regions ---
            genes_hit = ann_df['gene_id'].dropna().unique().tolist()

            # --- Debug: show number and few example genes ---
            print(f"🧠 Total unique genes hit: {len(genes_hit)}")
            if len(genes_hit) > 0:
                print("🔍 Example genes:", genes_hit[:10])
            else:
                print("⚠️ No genes found overlapping DMRs – check coordinate consistency or annotation file!")

            if len(genes_hit) > 0:
                # genes_hit is already a list of unique gene IDs
                enr = gp.profile(organism="athaliana", query=genes_hit)  # ⚠ change organism if needed
                enr.to_csv(f"{prefix}_pathway_enrichment.csv", index=False)
            else:
                print("⚠️ No genes found for enrichment. Skipping pathway analysis.")


            summary_list.append({
                "comparison": f"{a}_vs_{b}",
                "n_significant": len(diff_sig),
                "n_genes": len(genes_hit)
            })

    # Summary barplot
    summary_df = pd.DataFrame(summary_list)
    if not summary_df.empty:
        plt.figure(figsize=(5, 3))
        sns.barplot(data=summary_df, x="comparison", y="n_significant", color="steelblue")
        plt.title(f"Significant CpGs ({region})")
        plt.xticks(rotation=30)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{region}_summary_barplot.png"), dpi=300)
        plt.close()
    else:
        print(f"⚠ No significant comparisons found for {region}, skipping plot.")

    
    logging.info(f"Annotation dataframe shape: {ann_df.shape}")
    logging.info(f"Annotation columns: {list(ann_df.columns)}")
    logging.info(f"Example rows:\n{ann_df.head(5)}")
    logging.info(f"Unique genes hit: {len(genes_hit)}")
    logging.info(f"Example genes: {genes_hit[:10]}")
