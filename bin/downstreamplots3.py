import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from matplotlib_venn import venn2
import gzip
import numpy as np

# ===========================
# CONFIGURATION
# ===========================
input_dir = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/Artemia_Analysis_251028"
out_dir = "figures_output_AU"
os.makedirs(out_dir, exist_ok=True)

#barcode_pairs = [("21", "22"), ("23", "24")]
barcode_pairs = [("21", "24")]
regions = ["mitochondrial_AU", "nuclear"]
all_barcodes = ["21", "22", "23", "24"]

# ===========================
# HELPER FUNCTIONS
# ===========================
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
    wrapped = textwrap.fill(text, width=95)
    fig.text(0.5, -0.12, wrapped, ha='center', va='top', fontsize=9, wrap=True)

def safe_ls(pattern):
    return [p for p in os.popen(f"ls {pattern} 2>/dev/null").read().split() if os.path.exists(p)]

# ===========================
# MAIN ANALYSIS
# ===========================
captions = []

for region in regions:
    for b1, b2 in barcode_pairs:
        # --- Load methylation ---
        bed1 = safe_ls(f"{input_dir}/methylation/{region}/*barcode{b1}_mods.bed")
        bed2 = safe_ls(f"{input_dir}/methylation/{region}/*barcode{b2}_mods.bed")
        if not bed1 or not bed2:
            continue
        df1 = load_methylation_bed(bed1[0])
        df2 = load_methylation_bed(bed2[0])
        merged = pd.merge(df1, df2, on=["chrom", "start"], suffixes=(f"_{b1}", f"_{b2}"))
        merged["diff"] = merged[f"mod_prob_{b2}"] - merged[f"mod_prob_{b1}"]

        # === A. Methylation difference histogram ===
        plt.figure(figsize=(10, 4))
        sns.histplot(merged["diff"], bins=50, kde=True, color="teal")
        plt.title(f"Methylation Difference Distribution ({region} {b1}→{b2})")
        plt.xlabel("Δ Methylation Probability")
        plt.ylabel("Frequency")
        fig = plt.gcf()
        wrapped_caption(fig, f"Methylation difference between generations {b1} and {b2} for {region}.")
        fA = f"{out_dir}/{region}_meth_diff_{b1}_vs_{b2}.png"
        fig.savefig(fA, bbox_inches="tight", dpi=300)
        plt.close()

        # === B. Mutation overlap (Venn) ===
        vcf1 = safe_ls(f"{input_dir}/variants/nuclear/*barcode{b1}.vcf.gz")
        vcf2 = safe_ls(f"{input_dir}/variants/mitochondrial/*barcode{b2}.vcf.gz")
        if not vcf1 or not vcf2:
            continue
        mut1 = load_vcf_positions(vcf1[0])
        mut2 = load_vcf_positions(vcf2[0])
        set1 = set(zip(mut1["chrom"], mut1["pos"]))
        set2 = set(zip(mut2["chrom"], mut2["pos"]))

        plt.figure(figsize=(5, 5))
        venn2([set1, set2], set_labels=(f"{b1}", f"{b2}"))
        plt.title(f"Mutation Overlap ({region} {b1} vs {b2})")
        fig = plt.gcf()
        wrapped_caption(fig, f"Overlap of mutation loci between generations {b1} and {b2} in {region}.")
        fB = f"{out_dir}/{region}_mutation_venn_{b1}_vs_{b2}.png"
        fig.savefig(fB, bbox_inches="tight", dpi=300)
        plt.close()

        # === C. Combined methylation–mutation scatter ===
        mut_positions = pd.DataFrame(list(set1.union(set2)), columns=["chrom", "pos"])
        merged["pos"] = merged["start"]
        merged["chrom"] = merged["chrom"].astype(str)
        mut_positions["chrom"] = mut_positions["chrom"].astype(str)
        combined = pd.merge(merged, mut_positions, how="left", on=["chrom", "pos"])

        combined["mutation_present"] = combined["pos"].isin(mut_positions["pos"])

        plt.figure(figsize=(10, 4))
        sns.scatterplot(data=combined, x="pos", y="diff", hue="mutation_present", s=10, alpha=0.6)
        plt.title(f"Methylation vs Mutation Map ({region} {b1}→{b2})")
        plt.xlabel("Genomic Position")
        plt.ylabel("Δ Methylation")
        plt.legend(title="Mutation Present", loc="upper right")
        fig = plt.gcf()
        wrapped_caption(fig, f"Positions with methylation change and overlapping mutation in {region}.")
        fC = f"{out_dir}/{region}_meth_mut_scatter_{b1}_vs_{b2}.png"
        fig.savefig(fC, bbox_inches="tight", dpi=300)
        plt.close()

        # === D. Methylation → Mutation conversion ===
        overlap_sites = combined.query("mutation_present and abs(diff) > 0.2")
        plt.figure(figsize=(6, 4))
        sns.countplot(x=(overlap_sites['diff'] > 0).map({True: 'Gain', False: 'Loss'}))
        plt.title(f"Methylation→Mutation Conversions ({region} {b1}→{b2})")
        plt.xlabel("Methylation Change Direction")
        plt.ylabel("Count of Overlapping Mutations")
        fig = plt.gcf()
        wrapped_caption(fig, f"Counts of loci showing methylation change (>0.2) and mutation in {region}.")
        fD = f"{out_dir}/{region}_meth_to_mut_{b1}_vs_{b2}.png"
        fig.savefig(fD, bbox_inches="tight", dpi=300)
        plt.close()

        captions.append(f"{region.upper()} | {b1}→{b2}\nA: {fA}\nB: {fB}\nC: {fC}\nD: {fD}\n")

# === E. Trajectory 21→22→23→24 ===
traj_df = []
for region in regions:
    region_avgs = []
    for b in all_barcodes:
        bed = safe_ls(f"{input_dir}/methylation/{region}/*barcode{b}_mods.bed")
        if not bed:
            continue
        df = load_methylation_bed(bed[0])
        mean_meth = df["mod_prob"].mean()
        region_avgs.append({"region": region, "barcode": b, "avg_methylation": mean_meth})
    traj_df.extend(region_avgs)

traj_df = pd.DataFrame(traj_df)
print("\n[DEBUG] traj_df columns:", traj_df.columns.tolist())
print("[DEBUG] traj_df head:\n", traj_df.head())
plt.figure(figsize=(8, 5))
sns.lineplot(data=traj_df, x="barcode", y="avg_methylation", hue="region", marker="o")
plt.title("Methylation Trajectory (21→24)")
plt.xlabel("Generation")
plt.ylabel("Average Methylation Probability")
fig = plt.gcf()
wrapped_caption(fig, "Average methylation probability trend across all four generations.")
fE = f"{out_dir}/methylation_trajectory.png"
fig.savefig(fE, bbox_inches="tight", dpi=300)
plt.close()

# === F. Pathway summary (placeholder for GO) ===
pathways = pd.DataFrame({
    "Pathway": ["Oxidative phosphorylation", "ATP synthesis", "DNA repair", "Ion transport", "Stress response"],
    "Score": [0.95, 0.89, 0.76, 0.68, 0.61]
})
plt.figure(figsize=(8, 4))
sns.barplot(data=pathways, x="Score", y="Pathway", color="lightblue")
plt.title("Top Differentially Affected Pathways (Example)")
fig = plt.gcf()
wrapped_caption(fig, "Representative pathways potentially affected by methylation/mutation patterns.")
fF = f"{out_dir}/pathway_summary.png"
fig.savefig(fF, bbox_inches="tight", dpi=300)
plt.close()

# === G. Captions file ===
with open(f"{out_dir}/figure_captions.txt", "w") as f:
    f.write("\n\n".join(captions))
    f.write(f"\n\nE: {fE}\nF: {fF}\n")

print("✅ All figures A–G and captions generated successfully!")
