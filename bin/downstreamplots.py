#!/usr/bin/env python3
"""
plot_methylation_mutation.py

Usage:
  python plot_methylation_mutation.py --results <RESULTS_DIR> [--nuc-gff <NUC_GFF>] [--threads N]

Example:
  python plot_methylation_mutation.py --results /path/to/results --nuc-gff /path/to/assembly.gff

Outputs:
  results/figures/*.png
  results/figures/figure_captions.txt

Notes:
 - Script expects methylation files under:
     <RESULTS>/methylation/nuclear/*.bed
     <RESULTS>/methylation/mitochondrial_AP/*.bed
     <RESULTS>/methylation/mitochondrial_AU/*.bed
   and variant files under:
     <RESULTS>/variants/nuclear/*.vcf.gz
     <RESULTS>/variants/mitochondrial/*.vcf.gz  (if your pipeline put mito vcfs elsewhere, point RESULTS to parent)
 - The script auto-detects mod fraction column in modkit bed files (common positions: 5 or 11). If wrong, adjust `detect_mod_fraction_col()` logic below.
"""
import argparse
import os
import sys
from glob import glob
import gzip
import math
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from matplotlib_venn import venn2
from sklearn.preprocessing import scale
from scipy.cluster.hierarchy import linkage, leaves_list
plt.rcParams["figure.dpi"] = 150
sns.set(style="whitegrid")

BARCODES = {"21": "barcode21", "22": "barcode22", "23": "barcode23", "24": "barcode24"}
GEN_NAMES = {"21":"G1_1996","22":"G2_2001","23":"G3_2003","24":"G4_2005"}

def detect_mod_fraction_col(df):
    """
    Heuristic to find the mod fraction column in a modkit .bed file.
    Looks for a column with numeric values in [0,1] with decimals (e.g. 0.00..1.00).
    Returns 0-based column index.
    """
    for col in df.columns:
        try:
            vals = pd.to_numeric(df[col], errors="coerce")
            if vals.notna().sum() > 0:
                v = vals.dropna()
                if (v >= 0).all() and (v <= 1).all():
                    # a good candidate; prefer columns with decimals and >1 unique
                    if len(v.unique()) > 1:
                        return col
        except Exception:
            continue
    # fallback: use 4th (0-based 4) or 10th (0-based 10) if present
    if 4 in df.columns:
        return 4
    if 10 in df.columns:
        return 10
    raise ValueError("Could not detect mod fraction column automatically.")

def read_modkit_bed(path):
    """
    Read modkit .bed file into a DataFrame with columns:
    chr, start, end, extra..., mod_frac (float, between 0 and 1)
    """
    # read first 50 rows to detect columns
    df_head = pd.read_csv(path, sep="\t", header=None, nrows=200, engine="python")
    # full read (no headers), allow variable columns
    df = pd.read_csv(path, sep="\t", header=None, engine="python")
    # ensure at least 3 columns exist
    if df.shape[1] < 3:
        raise ValueError(f"{path} seems not to be a BED-like file (cols < 3).")
    # detect mod fraction column index
    col_idx = detect_mod_fraction_col(df)
    mod_frac = pd.to_numeric(df[col_idx], errors="coerce")
    out = pd.DataFrame({
        "chr": df[0].astype(str),
        "start": df[1].astype(int),
        "end": df[2].astype(int),
        "mod_frac": mod_frac
    })
    # drop NA mod_frac rows (if any)
    out = out.dropna(subset=["mod_frac"]).reset_index(drop=True)
    # create position key
    out["pos_key"] = out["chr"] + ":" + out["start"].astype(str)
    return out

def read_vcf_positions(vcf_path):
    """
    Return a set of position keys 'chr:start' for variants in the VCF (any ALT).
    """
    posset = set()
    try:
        vcf = pysam.VariantFile(vcf_path)
    except Exception as e:
        print(f"[WARN] Could not open VCF {vcf_path}: {e}", file=sys.stderr)
        return posset
    for rec in vcf.fetch():
        # pysam uses 1-based POS; our modkit start is 0- or 1-based? modkit bed start seems 481 for 481-482,
        # It looks 0-based start; vcf POS is 1-based. We'll compare using start (modkit start) +1 == VCF POS.
        key = f"{rec.chrom}:{rec.pos}"
        posset.add(key)
    return posset

def gather_sample_methylation(results_dir, category, barcode):
    """
    Search results/methylation/<category> for files containing barcode string,
    read the modkit bed and return DataFrame.
    """
    d = os.path.join(results_dir, "methylation", category)
    patterns = [os.path.join(d, f"*{barcode}*.bed"), os.path.join(d, f"*{barcode}*_mods.bed")]
    files = []
    for p in patterns:
        files.extend(sorted(glob(p)))
    if not files:
        return None, None
    # choose first match (should be unique)
    path = files[0]
    df = read_modkit_bed(path)
    return df, path

def gather_sample_vcf(results_dir, category, barcode):
    """
    Try typical vcf locations. category: 'nuclear' or 'mitochondrial' (we search variants/nuclear and variants/mitochondrial).
    """
    variants_dir = os.path.join(results_dir, "variants")
    if category == "nuclear":
        d = os.path.join(variants_dir, "nuclear")
    else:
        # many pipelines put mito vcfs under variants/mitochondrial
        d = os.path.join(variants_dir, "mitochondrial")
    patterns = [os.path.join(d, f"*{barcode}*.vcf.gz"), os.path.join(d, f"*{barcode}*.vcf")]
    files = []
    for p in patterns:
        files.extend(sorted(glob(p)))
    if not files:
        return None
    return files[0]

def scatter_methylation(dfA, dfB, labelA, labelB, outpng):
    merged = pd.merge(dfA, dfB, on="pos_key", suffixes=(f"_{labelA}", f"_{labelB}"))
    if merged.empty:
        print(f"[WARN] No overlapping CpGs between {labelA} and {labelB} to scatter.")
        return merged
    plt.figure(figsize=(6,6))
    sns.scatterplot(x=f"mod_frac_{labelA}", y=f"mod_frac_{labelB}", data=merged, s=10, alpha=0.5)
    plt.plot([0,1],[0,1], '--', color='red')
    plt.xlabel(f"{GEN_NAMES[labelA]} methylation")
    plt.ylabel(f"{GEN_NAMES[labelB]} methylation")
    plt.title(f"Methylation: {GEN_NAMES[labelA]} vs {GEN_NAMES[labelB]}")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()
    return merged

def histogram_delta(merged, labelA, labelB, outpng):
    merged["delta"] = merged[f"mod_frac_{labelB}"] - merged[f"mod_frac_{labelA}"]
    plt.figure(figsize=(6,4))
    sns.histplot(merged["delta"], bins=50, kde=False)
    plt.axvline(0, color='red', linestyle='--')
    plt.title(f"Δ methylation {GEN_NAMES[labelB]} - {GEN_NAMES[labelA]}")
    plt.xlabel("Δ methylation fraction")
    plt.ylabel("CpG count")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def venn_methylation(dfA, dfB, labelA, labelB, outpng):
    setA = set(dfA["pos_key"])
    setB = set(dfB["pos_key"])
    plt.figure(figsize=(5,5))
    venn2([setA, setB], set_labels=(GEN_NAMES[labelA], GEN_NAMES[labelB]))
    plt.title(f"Venn: methylation CpGs {GEN_NAMES[labelA]} vs {GEN_NAMES[labelB]}")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()
    return setA, setB

def plot_mutation_comparison(vcfA, vcfB, labelA, labelB, outpng):
    setA = read_vcf_positions(vcfA) if vcfA else set()
    setB = read_vcf_positions(vcfB) if vcfB else set()
    shared = setA.intersection(setB)
    uniqueA = len(setA - setB)
    uniqueB = len(setB - setA)
    sharedN = len(shared)
    labels = [f"Unique {GEN_NAMES[labelA]}", "Shared", f"Unique {GEN_NAMES[labelB]}"]
    vals = [uniqueA, sharedN, uniqueB]
    plt.figure(figsize=(5,4))
    sns.barplot(x=labels, y=vals)
    plt.ylabel("Variant count")
    plt.title(f"Mutation comparison {GEN_NAMES[labelA]} vs {GEN_NAMES[labelB]}")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()
    return setA, setB

def methylation_to_mutation_overlap(df_from, vcf_to, label_from, label_to, outpng):
    """
    Are methylated CpGs in generation 'from' overlapping variant positions in generation 'to'?
    """
    if df_from is None or vcf_to is None:
        print(f"[INFO] Missing data for methyl->mut check {label_from}->{label_to}")
        return 0,0
    meth_positions = set(df_from["pos_key"])
    mut_positions = read_vcf_positions(vcf_to)
    # transform meth pos (0-based start) to 1-based pos for matching VCF POS:
    # If modkit start is 0-based start, then VCF pos equals start+1. We will try both possibilities:
    meth_positions_v1 = set()
    for k in meth_positions:
        chrom, s = k.split(":")
        try:
            pos0 = int(s)
        except:
            continue
        # two possible keys in vcf sets are chrom:pos where pos is 1-based
        meth_positions_v1.add(f"{chrom}:{pos0+1}")  # assume modkit 0-based
        meth_positions_v1.add(f"{chrom}:{pos0}")    # in case modkit used 1-based
    overlap = meth_positions_v1.intersection(mut_positions)
    n_overlap = len(overlap)
    # Plot simple counts
    plt.figure(figsize=(5,3))
    sns.barplot(x=["methylated in "+GEN_NAMES[label_from], "mutated in "+GEN_NAMES[label_to], "overlap"],
                y=[len(meth_positions), len(mut_positions), n_overlap])
    plt.xticks(rotation=45, ha='right')
    plt.title(f"Methylation( {GEN_NAMES[label_from]} ) -> Mutation( {GEN_NAMES[label_to]} ) overlap")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()
    return len(meth_positions), n_overlap

def heatmap_top_cpgs(all_df_dict, outpng, top_n=500):
    """
    all_df_dict: dict barcode->df with pos_key and mod_frac
    Build matrix for union of top variable CpGs and cluster.
    """
    # build DataFrame with rows=pos_key, cols=sample
    merged = None
    for b, df in all_df_dict.items():
        if df is None:
            continue
        tmp = df[["pos_key", "mod_frac"]].set_index("pos_key").rename(columns={"mod_frac": b})
        if merged is None:
            merged = tmp
        else:
            merged = merged.join(tmp, how="outer")
    if merged is None or merged.shape[0] == 0:
        print("[WARN] No methylation data found for heatmap.")
        return
    merged = merged.fillna(0)
    # choose top variable rows
    merged['var'] = merged.var(axis=1)
    top = merged.sort_values("var", ascending=False).head(top_n).drop(columns="var")
    # cluster rows
    Z = linkage(scale(top), method='average')
    row_order = leaves_list(Z)
    top_ord = top.iloc[row_order, :]
    plt.figure(figsize=(8, max(6, top_ord.shape[0]/50)))
    sns.heatmap(top_ord, cmap="viridis", cbar_kws={"label":"methylation fraction"})
    plt.title("Top variable CpGs across generations")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, help="Path to RESULTS directory (the one containing methylation/ and variants/)")
    parser.add_argument("--nuc-gff", default=None, help="Optional path to nuclear GFF for gene annotation (not required)")
    parser.add_argument("--outdir", default=None, help="Output directory (default: <RESULTS>/figures)")
    args = parser.parse_args()

    RESULTS = args.results
    OUTDIR = args.outdir or os.path.join(RESULTS, "figures")
    os.makedirs(OUTDIR, exist_ok=True)

    # containers for data
    data = {
        "mito_AP": {},
        "mito_AU": {},
        "nuclear": {}
    }
    vcf_paths = {
        "mito": {},
        "nuclear": {}
    }

    # load methylation tables for each barcode and each category
    for cat in ["mitochondrial_AP", "mitochondrial_AU", "nuclear"]:
        short_cat = "mito_AP" if "AP" in cat else ("mito_AU" if "AU" in cat else "nuclear")
        for bc in ["21","22","23","24"]:
            df, path = gather_sample_methylation(RESULTS, cat if cat!="nuclear" else "nuclear", BARCODES[bc])
            data[short_cat][bc] = df
    # load vcf paths
    for bc in ["21","22","23","24"]:
        v21 = gather_sample_vcf(RESULTS, "nuclear", BARCODES[bc])
        v_mito = gather_sample_vcf(RESULTS, "mitochondrial", BARCODES[bc])
        vcf_paths["nuclear"][bc] = v21
        vcf_paths["mito"][bc] = v_mito

    captions = []
    # 1) Mitochondrial comparisons: 21 vs 22 (AP) and 23 vs 24 (AU) -> methylation
    # -> 21 vs 22 (mito_AP)
    print("Plotting mitochondrial comparisons 21 vs 22 (AP)...")
    df21 = data["mito_AP"].get("21")
    df22 = data["mito_AP"].get("22")
    if df21 is not None and df22 is not None:
        m_merged = scatter_methylation(df21, df22, "21", "22", os.path.join(OUTDIR, "mito_AP_21_vs_22_scatter.png"))
        if isinstance(m_merged, pd.DataFrame) and not m_merged.empty:
            histogram_delta(m_merged, "21", "22", os.path.join(OUTDIR, "mito_AP_21_vs_22_delta_hist.png"))
        venn_methylation(df21, df22, "21", "22", os.path.join(OUTDIR, "mito_AP_21_vs_22_venn.png"))
        captions.append(("Fig: mito_AP_21_vs_22_scatter.png",
                         "Comparison of mitochondrial CpG methylation between generation 21 (1996) and 22 (2001). Scatter of methylation fractions per CpG."))
    else:
        print("[WARN] Missing mito_AP data for 21 or 22.")

    # 2) mito AU 23 vs 24
    print("Plotting mitochondrial comparisons 23 vs 24 (AU)...")
    df23 = data["mito_AU"].get("23")
    df24 = data["mito_AU"].get("24")
    if df23 is not None and df24 is not None:
        m_merged = scatter_methylation(df23, df24, "23", "24", os.path.join(OUTDIR, "mito_AU_23_vs_24_scatter.png"))
        if isinstance(m_merged, pd.DataFrame) and not m_merged.empty:
            histogram_delta(m_merged, "23", "24", os.path.join(OUTDIR, "mito_AU_23_vs_24_delta_hist.png"))
        venn_methylation(df23, df24, "23", "24", os.path.join(OUTDIR, "mito_AU_23_vs_24_venn.png"))
        captions.append(("Fig: mito_AU_23_vs_24_scatter.png",
                         "Comparison of mitochondrial CpG methylation between generation 23 (2003) and 24 (2005). Scatter of methylation fractions per CpG."))
    else:
        print("[WARN] Missing mito_AU data for 23 or 24.")

    # 3) Mutations: mito 21 vs 22 and 23 vs 24
    print("Plotting mitochondrial mutation comparisons...")
    v21 = vcf_paths["mito"].get("21")
    v22 = vcf_paths["mito"].get("22")
    if v21 or v22:
        plot_mutation_comparison(v21, v22, "21", "22", os.path.join(OUTDIR, "mito_21_vs_22_mutation_counts.png"))
        captions.append(("Fig: mito_21_vs_22_mutation_counts.png",
                         "Mitochondrial variant comparison between generation 21 and 22. Bars show unique and shared variant counts."))
    v23 = vcf_paths["mito"].get("23")
    v24 = vcf_paths["mito"].get("24")
    if v23 or v24:
        plot_mutation_comparison(v23, v24, "23", "24", os.path.join(OUTDIR, "mito_23_vs_24_mutation_counts.png"))
        captions.append(("Fig: mito_23_vs_24_mutation_counts.png",
                         "Mitochondrial variant comparison between generation 23 and 24. Bars show unique and shared variant counts."))

    # 4) Nuclear comparisons: methylation and mutation for 21 vs 22, 23 vs 24
    print("Plotting nuclear methylation and mutation comparisons...")
    dn21 = data["nuclear"].get("21")
    dn22 = data["nuclear"].get("22")
    dn23 = data["nuclear"].get("23")
    dn24 = data["nuclear"].get("24")
    if dn21 is not None and dn22 is not None:
        merged = scatter_methylation(dn21, dn22, "21", "22", os.path.join(OUTDIR,"nuclear_21_vs_22_scatter.png"))
        if isinstance(merged, pd.DataFrame) and not merged.empty:
            histogram_delta(merged, "21", "22", os.path.join(OUTDIR,"nuclear_21_vs_22_delta_hist.png"))
        venn_methylation(dn21, dn22, "21", "22", os.path.join(OUTDIR,"nuclear_21_vs_22_venn.png"))
        captions.append(("Fig: nuclear_21_vs_22_scatter.png", "Nuclear CpG methylation: Generation 21 vs 22."))
    else:
        print("[WARN] Missing nuclear methylation for 21 or 22.")

    if dn23 is not None and dn24 is not None:
        merged = scatter_methylation(dn23, dn24, "23", "24", os.path.join(OUTDIR,"nuclear_23_vs_24_scatter.png"))
        if isinstance(merged, pd.DataFrame) and not merged.empty:
            histogram_delta(merged, "23", "24", os.path.join(OUTDIR,"nuclear_23_vs_24_delta_hist.png"))
        venn_methylation(dn23, dn24, "23", "24", os.path.join(OUTDIR,"nuclear_23_vs_24_venn.png"))
        captions.append(("Fig: nuclear_23_vs_24_scatter.png", "Nuclear CpG methylation: Generation 23 vs 24."))
    else:
        print("[WARN] Missing nuclear methylation for 23 or 24.")

    # nuclear mutation comparisons
    print("Plotting nuclear mutation comparisons...")
    nv21 = vcf_paths["nuclear"].get("21")
    nv22 = vcf_paths["nuclear"].get("22")
    nv23 = vcf_paths["nuclear"].get("23")
    nv24 = vcf_paths["nuclear"].get("24")
    if nv21 or nv22:
        plot_mutation_comparison(nv21, nv22, "21", "22", os.path.join(OUTDIR,"nuclear_21_vs_22_mutation_counts.png"))
        captions.append(("Fig: nuclear_21_vs_22_mutation_counts.png","Nuclear variant comparison 21 vs 22."))
    if nv23 or nv24:
        plot_mutation_comparison(nv23, nv24, "23", "24", os.path.join(OUTDIR,"nuclear_23_vs_24_mutation_counts.png"))
        captions.append(("Fig: nuclear_23_vs_24_mutation_counts.png","Nuclear variant comparison 23 vs 24."))

    # 5) Methylation -> Mutation checks
    print("Checking methylation -> mutation overlaps...")
    mm21_22 = methylation_to_mutation_overlap(dn21, nv22, "21", "22", os.path.join(OUTDIR,"nuc_meth21_to_mut22_overlap.png"))
    mm23_24 = methylation_to_mutation_overlap(dn23, nv24, "23", "24", os.path.join(OUTDIR,"nuc_meth23_to_mut24_overlap.png"))
    captions.append(("Fig: nuc_meth21_to_mut22_overlap.png","Nuclear methylation (21) → mutation (22) overlap counts."))
    captions.append(("Fig: nuc_meth23_to_mut24_overlap.png","Nuclear methylation (23) → mutation (24) overlap counts."))

    # mito methylation -> mutation
    mmm21_22 = methylation_to_mutation_overlap(df21, v22, "21", "22", os.path.join(OUTDIR,"mito_meth21_to_mut22_overlap.png")) if (df21 is not None) else (0,0)
    mmm23_24 = methylation_to_mutation_overlap(df23, v24, "23", "24", os.path.join(OUTDIR,"mito_meth23_to_mut24_overlap.png")) if (df23 is not None) else (0,0)
    captions.append(("Fig: mito_meth21_to_mut22_overlap.png","Mitochondrial methylation (21) → mutation (22) overlap counts."))
    captions.append(("Fig: mito_meth23_to_mut24_overlap.png","Mitochondrial methylation (23) → mutation (24) overlap counts."))

    # 6) Multi-generation heatmap (top variable CpGs)
    print("Building heatmap across generations (mitochondrial and nuclear separately)...")
    # mitochondrial combined union from AP and AU: use whichever exists
    mito_all = {}
    for bc in ["21","22","23","24"]:
        # prefer mito_AP for 21/22 and mito_AU for 23/24
        if bc in ["21","22"]:
            mito_all[bc] = data["mito_AP"].get(bc) or data["mito_AU"].get(bc)
        else:
            mito_all[bc] = data["mito_AU"].get(bc) or data["mito_AP"].get(bc)
    heatmap_top_cpgs(mito_all, os.path.join(OUTDIR,"mito_heatmap_top_variable.png"))
    captions.append(("Fig: mito_heatmap_top_variable.png","Heatmap of top variable mitochondrial CpGs across generations."))

    heatmap_top_cpgs(data["nuclear"], os.path.join(OUTDIR,"nuclear_heatmap_top_variable.png"))
    captions.append(("Fig: nuclear_heatmap_top_variable.png","Heatmap of top variable nuclear CpGs across generations."))

    # 7) write captions file
    with open(os.path.join(OUTDIR,"figure_captions.txt"), "w") as fh:
        fh.write("Figure captions generated by plot_methylation_mutation.py\n\n")
        for name, txt in captions:
            fh.write(f"{name}\n{txt}\n\n")
    print("All figures/written to", OUTDIR)
    print("Finished.")

if __name__ == "__main__":
    main()
