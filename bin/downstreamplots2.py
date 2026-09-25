import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from matplotlib_venn import venn2
import gzip

# --- CONFIG ---
input_dir = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results"
out_dir = "figures_output"
os.makedirs(out_dir, exist_ok=True)

barcode_pairs = [("21", "22"), ("23", "24")]
regions = ["mitochondrial_AP", "nuclear"]

# --- Helper functions ---
def load_methylation_bed(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "mod", "score", "strand",
                            "col7", "col8", "color", "col10", "mod_prob", 
                            "col12", "col13", "col14", "col15", "col16", "col17"])
    df = df[["chrom", "start", "mod", "mod_prob"]]
    return df

def load_vcf_positions(vcf_path):
    positions = []
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, *_ = line.split("\t")
            positions.append((chrom, int(pos)))
    return pd.DataFrame(positions, columns=["chrom", "pos"])

def wrapped_caption(fig, text):
    wrapped = textwrap.fill(text, width=90)
    fig.text(0.5, -0.1, wrapped, ha='center', va='top', fontsize=9, wrap=True)

# --- Analysis loop ---
captions = []

for region in regions:
    for b1, b2 in barcode_pairs:
        # --- Load methylation data ---
        bed1 = f"{input_dir}/methylation/{region}/*barcode{b1}_mods.bed"
        bed2 = f"{input_dir}/methylation/{region}/*barcode{b2}_mods.bed"
        paths1 = [p for p in os.popen(f"ls {bed1}").read().split() if os.path.exists(p)]
        paths2 = [p for p in os.popen(f"ls {bed2}").read().split() if os.path.exists(p)]
        if not paths1 or not paths2:
            continue
        df1 = load_methylation_bed(paths1[0])
        df2 = load_methylation_bed(paths2[0])

        # --- Merge methylation ---
        merged = pd.merge(df1, df2, on=["chrom", "start"], suffixes=(f"_{b1}", f"_{b2}"))
        merged["diff"] = merged[f"mod_prob_{b2}"] - merged[f"mod_prob_{b1}"]

        # --- Heatmap of methylation difference ---
        plt.figure(figsize=(10, 4))
        sns.histplot(merged["diff"], bins=50, kde=True, color="teal")
        plt.title(f"Methylation Difference Distribution ({region} {b1}→{b2})")
        plt.xlabel("Δ Methylation Probability")
        plt.ylabel("Frequency")

        fig = plt.gcf()
        wrapped_caption(
            fig,
            f"Figure: Methylation difference between generations {b1} and {b2} "
            f"for {region}. Positive Δ indicates higher methylation in {b2}."
        )

        fig.savefig(f"{out_dir}/{region}_meth_diff_{b1}_vs_{b2}.png", bbox_inches="tight", dpi=300)
        plt.close()

        # --- Mutation comparison ---
        vcf1 = f"{input_dir}/variants/nuclear/*barcode{b1}.vcf.gz"
        vcf2 = f"{input_dir}/variants/mitochondrial/*barcode{b2}.vcf.gz"
        paths1 = [p for p in os.popen(f"ls {vcf1}").read().split() if os.path.exists(p)]
        paths2 = [p for p in os.popen(f"ls {vcf2}").read().split() if os.path.exists(p)]
        if not paths1 or not paths2:
            continue
        mut1 = load_vcf_positions(paths1[0])
        mut2 = load_vcf_positions(paths2[0])

        set1 = set(zip(mut1["chrom"], mut1["pos"]))
        set2 = set(zip(mut2["chrom"], mut2["pos"]))

        # --- Venn diagram ---
        plt.figure(figsize=(5, 5))
        venn2([set1, set2], set_labels=(f"{b1}", f"{b2}"))
        plt.title(f"Mutation Overlap ({region} {b1} vs {b2})")

        fig = plt.gcf()
        wrapped_caption(
            fig,
            f"Figure: Overlap of mutation positions between generations {b1} and {b2} "
            f"for {region}. Unique and shared variant loci are shown."
        )

        fig.savefig(f"{out_dir}/{region}_mutation_venn_{b1}_vs_{b2}.png", bbox_inches="tight", dpi=300)
        plt.close()

        # --- Save caption summary ---
        captions.append(
            f"{region.upper()} | {b1} vs {b2}\n"
            f"1. Methylation difference: {out_dir}/{region}_meth_diff_{b1}_vs_{b2}.png\n"
            f"2. Mutation overlap: {out_dir}/{region}_mutation_venn_{b1}_vs_{b2}.png\n"
        )

# --- Write caption file ---
with open(f"{out_dir}/figure_captions.txt", "w") as f:
    f.write("\n\n".join(captions))

print("✅ All figures and captions generated successfully!")
