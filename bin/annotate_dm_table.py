import os
import pandas as pd
from glob import glob

# --- Paths ---
base_dir = "/data/Jyotirmoy/CF-Projects/External/Parisa/CFFMHS-QU-12/artemia/artemia_analysis_251021/results"
regions = ["nuclear", "mitochondrial_AP", "mitochondrial_AU"]

# --- File size threshold in bytes (50 KB) ---
MIN_SIZE = 50 * 1024

for region in regions:
    print(f"\n🧬 Processing region: {region}")

    methyl_dir = os.path.join(base_dir, "methylation", region)
    diff_dir = os.path.join(base_dir, "differential_pathway_results", region)
    os.makedirs(diff_dir, exist_ok=True)

    bed_files = glob(os.path.join(methyl_dir, "*_mods.bed"))
    if not bed_files:
        print(f"⚠️ No BED files found in {methyl_dir}")
        continue

    coord_list = []
    skipped = 0
    for bed in bed_files:
        size = os.path.getsize(bed)
        if size <= MIN_SIZE:
            print(f"⚠️ Skipping small/non-standard file: {os.path.basename(bed)} ({size/1024:.1f} KB)")
            skipped += 1
            continue

        try:
            df = pd.read_csv(bed, sep="\t", header=None)
            if df.empty:
                print(f"⚠️ Skipping empty (parsed) file: {os.path.basename(bed)}")
                skipped += 1
                continue
            # keep first 3 columns
            df = df.iloc[:, :3]
            df.columns = ["chrom", "start", "end"]
            df["pos"] = df["start"].astype(int)
            coord_list.append(df[["chrom", "pos"]])
        except Exception as e:
            print(f"❌ Failed to read {os.path.basename(bed)}: {e}")
            skipped += 1

    if not coord_list:
        print(f"⚠️ No valid coordinate data for region {region} (all files skipped)")
        continue

    coord_df = pd.concat(coord_list).drop_duplicates()
    print(f"✅ Collected {len(coord_df)} unique positions from {len(bed_files) - skipped}/{len(bed_files)} valid BED files")

    # --- Differential tables ---
    diff_files = glob(os.path.join(diff_dir, "*_differential.csv"))
    if not diff_files:
        print(f"⚠️ No differential files found in {diff_dir}")
        continue

    for diff_file in diff_files:
        name = os.path.basename(diff_file)
        print(f"🔍 Annotating {name}")

        diff_df = pd.read_csv(diff_file)
        if "pos" not in diff_df.columns:
            diff_df.reset_index(inplace=True)
            diff_df.rename(columns={"index": "pos"}, inplace=True)

        merged = pd.merge(diff_df, coord_df, on="pos", how="left")
        out_path = os.path.join(diff_dir, name.replace(".csv", "_annotated.csv"))
        merged.to_csv(out_path, index=False)
        print(f"💾 Saved: {out_path}")
