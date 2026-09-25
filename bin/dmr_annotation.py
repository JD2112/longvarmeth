#!/usr/bin/env python3
"""
Annotate differential methylation CSVs with GFF gene overlaps.

- Finds files in: differential_pathway_results/*_differential.csv
- Uses:
    nuclear GFF: GCF_032884065.1_ASM3288406v1_genomic.gff
    mitochondrial GFF: MN240408_1.gff3
- Decides which to use based on filename (contains "mitochondrial" → mitochondrial GFF)
- Outputs: *_differential_annotated.csv next to each input CSV.

Requires: pybedtools (recommended) OR pure Python fallback.
"""

import os
import sys
from glob import glob
import pandas as pd

# ======================= CONFIG =======================
BASE_DIR = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia"
DIFF_GLOB = os.path.join(BASE_DIR, "artemia_analysis_251021/results/differential_pathway_results", "*_differential.csv")
NUCLEAR_GFF = os.path.join(BASE_DIR, "reference_genome", "GCF_032884065.1_ASM3288406v1_genomic.gff")
MITO_GFF = os.path.join(BASE_DIR, "reference_genome", "MN240408_1.gff3")
OUT_SUFFIX = "_annotated.csv"
TMPBED_DIR = os.path.join(BASE_DIR, "tmp_annotate_beds")

os.makedirs(TMPBED_DIR, exist_ok=True)

# ======================= UTILS =======================
def find_diff_files(pattern=DIFF_GLOB):
    files = sorted(glob(pattern))
    if not files:
        print(f"[ERROR] No differential files found with pattern: {pattern}")
    return files

# try to import pybedtools
USE_PYBEDTOOLS = False
try:
    import pybedtools
    USE_PYBEDTOOLS = True
    print("[INFO] pybedtools found — using fast bed intersection.")
except Exception:
    print("[INFO] pybedtools NOT found — using pure-Python overlap (slower).")

def load_gff_genes(gff_path):
    """Parse GFF and return DataFrame with columns: chrom, start, end, gene_id, gene_name"""
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
            # only gene-like features
            if ftype.lower() not in ("gene", "mrna", "cds", "transcript", "exon"):
                continue
            gene_id = None
            gene_name = None
            for token in attrs.split(";"):
                token = token.strip()
                if token.startswith("ID=") and not gene_id:
                    gene_id = token.split("=", 1)[1]
                if token.startswith("Name=") and not gene_name:
                    gene_name = token.split("=", 1)[1]
                if token.lower().startswith("gene=") and not gene_name:
                    gene_name = token.split("=", 1)[1]
                if token.lower().startswith("gene_id=") and not gene_id:
                    gene_id = token.split("=", 1)[1]
            if not gene_id and gene_name:
                gene_id = gene_name
            if not gene_id:
                first = attrs.split(";")[0].strip()
                gene_id = first.split("=")[-1] if "=" in first else first
            genes.append((chrom, int(start), int(end), gene_id, gene_name or ""))
    gdf = pd.DataFrame(genes, columns=["chrom", "start", "end", "gene_id", "gene_name"])
    gdf = gdf.drop_duplicates(subset=["chrom", "start", "end", "gene_id"])
    print(f"[INFO] Loaded {len(gdf)} features from {os.path.basename(gff_path)}")
    return gdf

def ensure_chrom_pos(df):
    """Ensure 'chrom' and 'pos' columns exist and standardized"""
    dfc = df.copy()
    candidates = {
        "chrom": ["chrom", "chr", "Chrom", "CHROM"],
        "pos": ["pos", "start", "position", "POS", "Start"]
    }
    found = {}
    for key, opts in candidates.items():
        for o in opts:
            if o in dfc.columns:
                found[key] = o
                break
    if "chrom" not in found:
        raise ValueError("No chromosome column found. Expected one of chrom/chr/Chrom/CHROM.")
    if "pos" not in found:
        dfc = dfc.reset_index().rename(columns={"index": "pos"})
        found["pos"] = "pos"
    dfc = dfc.rename(columns={found["chrom"]: "chrom", found["pos"]: "pos"})
    return dfc

def annotate_with_pybedtools(diff_df, gff_path, basename):
    """Use pybedtools to intersect diff positions with GFF genes."""
    bed_path = os.path.join(TMPBED_DIR, f"{basename}.bed")
    df = diff_df.copy()
    df["pos"] = df["pos"].astype(int)
    df["bed_start"] = df["pos"]
    df["bed_end"] = df["pos"] + 1
    df[["chrom", "bed_start", "bed_end"]].to_csv(bed_path, sep="\t", header=False, index=False)

    bt_pos = pybedtools.BedTool(bed_path)
    bt_gff = pybedtools.BedTool(gff_path)
    inter = bt_pos.intersect(bt_gff, wa=True, wb=True)

    hits = {}
    for line in inter:
        a = list(line)
        pos_chrom = a[0]
        pos_start = int(a[1])
        gff_attrs = a[-1]
        gene_id = ""
        gene_name = ""
        for token in gff_attrs.split(";"):
            t = token.strip()
            if t.startswith("ID="):
                gene_id = t.split("=", 1)[1]
            if t.startswith("Name="):
                gene_name = t.split("=", 1)[1]
            if t.lower().startswith("gene=") and not gene_name:
                gene_name = t.split("=", 1)[1]
        if not gene_id and gene_name:
            gene_id = gene_name
        key = (pos_chrom, pos_start)
        hits.setdefault(key, []).append((gene_id, gene_name))

    gene_id_list, gene_name_list = [], []
    for _, row in df.iterrows():
        key = (str(row["chrom"]), int(row["bed_start"]))
        if key in hits:
            gids = [g[0] for g in hits[key] if g[0]]
            gnames = [g[1] for g in hits[key] if g[1]]
            gene_id_list.append(",".join(sorted(set(gids))) if gids else "")
            gene_name_list.append(",".join(sorted(set(gnames))) if gnames else "")
        else:
            gene_id_list.append("")
            gene_name_list.append("")
    df["gene_id"] = gene_id_list
    df["gene_name"] = gene_name_list
    return df.drop(columns=["bed_start", "bed_end"])

def annotate_with_pure_python(diff_df, gff_df):
    """Pure-Python overlap fallback."""
    df = diff_df.copy()
    df["pos"] = df["pos"].astype(int)
    df["gene_id"], df["gene_name"] = "", ""
    ggroup = {chrom: group.sort_values("start")[["start", "end", "gene_id", "gene_name"]].values.tolist()
              for chrom, group in gff_df.groupby("chrom")}
    for i, row in df.iterrows():
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        if chrom not in ggroup:
            continue
        hits = [(gid, gname) for start, end, gid, gname in ggroup[chrom] if start <= pos <= end]
        if hits:
            gids = [h[0] for h in hits if h[0]]
            gnames = [h[1] for h in hits if h[1]]
            df.at[i, "gene_id"] = ",".join(sorted(set(gids))) if gids else ""
            df.at[i, "gene_name"] = ",".join(sorted(set(gnames))) if gnames else ""
    return df

def annotate_file(diff_path, gff_df, gff_path):
    """Annotate one CSV file."""
    print(f"[INFO] Annotating: {diff_path} using {os.path.basename(gff_path)}")
    df = pd.read_csv(diff_path)
    df = ensure_chrom_pos(df)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").dropna().astype(int)
    base = os.path.basename(diff_path)
    name = os.path.splitext(base)[0]
    if USE_PYBEDTOOLS:
        ann = annotate_with_pybedtools(df, gff_path, name)
    else:
        ann = annotate_with_pure_python(df, gff_df)
    outpath = diff_path.replace("_differential.csv", "_differential_annotated.csv")
    ann.to_csv(outpath, index=False)
    print(f"[OK] Wrote annotated file: {outpath}")

# ======================= MAIN =======================
def main():
    diff_files = find_diff_files()
    if not diff_files:
        return

    # Pre-load both GFFs for efficiency
    print("[INFO] Loading nuclear and mitochondrial annotations...")
    gff_nuclear = load_gff_genes(NUCLEAR_GFF)
    gff_mito = load_gff_genes(MITO_GFF)

    for fp in diff_files:
        try:
            if "mitochondrial" in os.path.basename(fp).lower():
                annotate_file(fp, gff_mito, MITO_GFF)
            else:
                annotate_file(fp, gff_nuclear, NUCLEAR_GFF)
        except Exception as e:
            print(f"[ERROR] Failed to annotate {fp}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
