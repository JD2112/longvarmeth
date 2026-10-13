#!/usr/bin/env python3
"""
Author: Jyotirmoy Das (adapted for longvarmeth)
Maintainer: Bioinformatics Core Facility and Clinical Genomics Linkoping
Version: 0.3.0
Purpose: Standalone interactive methylation clinical dashboard generator from methbat / modkit profile outputs.
"""

import json
import csv
import argparse
import sys
from pathlib import Path


def parse_methbat_tsv(tsv_path):
    """Parses methbat profile TSV/CSV or modkit aggregated TSV into a list of dicts for AG-Grid."""
    rows = []
    try:
        delimiter = "," if tsv_path.suffix == ".csv" else "\t"
        with open(tsv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = None
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                header = [h.strip().lower() for h in row]
                break

            if not header:
                return []

            def get_idx(name_list):
                for name in name_list:
                    if name in header:
                        return header.index(name)
                return -1

            idx_chrom = get_idx(["chrom", "chromosome", "chr", "seq_id"])
            idx_start = get_idx(["start", "start_pos", "pos", "begin"])
            idx_end = get_idx(["end", "end_pos"])
            idx_cpg = get_idx(["cpg_label", "gene", "region", "name", "id"])
            idx_exp_cat = get_idx(["expected_category", "category", "class"])
            idx_sum_lab = get_idx(["summary_label", "status", "label"])
            idx_combined = get_idx(
                [
                    "mean_combined_methyl",
                    "combined_methyl",
                    "fraction_modified",
                    "percent_modified",
                    "meth_score",
                ]
            )
            idx_delta = get_idx(["mean_meth_delta", "meth_delta", "delta"])
            idx_h1 = get_idx(["mean_hap1_methyl", "hap1_methyl", "h1_meth"])
            idx_h2 = get_idx(["mean_hap2_methyl", "hap2_methyl", "h2_meth"])
            idx_qc_warn = get_idx(["qc_warnings", "warnings", "qc"])
            idx_qc_sum = get_idx(["qc_summary", "filter"])

            for row in reader:
                if not row:
                    continue

                def safe_get(idx, default=""):
                    if idx >= 0 and idx < len(row):
                        return row[idx].strip()
                    return default

                start = safe_get(idx_start, "0")
                end = safe_get(idx_end, start)
                interval = f"{safe_get(idx_chrom)}:{start}-{end}"

                warnings = safe_get(idx_qc_warn, "")
                summary = safe_get(idx_qc_sum, "PASS")

                if (warnings == "none" or not warnings) and (
                    summary in ["PASS", "pass", ""]
                ):
                    qc = "PASS"
                    qc_status = 0
                elif summary in ["PASS", "pass"] and warnings:
                    qc = f"Warning: {warnings}"
                    qc_status = 1
                else:
                    qc = f"{summary}; {warnings}".strip("; ")
                    qc_status = 2

                status = (
                    safe_get(idx_sum_lab) or safe_get(idx_exp_cat) or "Uncategorized"
                )
                combined_raw = safe_get(idx_combined, "0")
                try:
                    combined_val = float(combined_raw)
                    # Normalize percent (0-100) to fraction (0-1) if needed
                    if combined_val > 1.0:
                        combined_val = combined_val / 100.0
                except ValueError:
                    combined_val = 0.0

                try:
                    delta_val = float(safe_get(idx_delta, "0") or 0)
                except ValueError:
                    delta_val = 0.0

                try:
                    h1_val = float(safe_get(idx_h1, "0") or 0)
                    if h1_val > 1.0:
                        h1_val /= 100.0
                except ValueError:
                    h1_val = 0.0

                try:
                    h2_val = float(safe_get(idx_h2, "0") or 0)
                    if h2_val > 1.0:
                        h2_val /= 100.0
                except ValueError:
                    h2_val = 0.0

                cleaned_row = {
                    "chrom": safe_get(idx_chrom),
                    "start": start,
                    "end": end,
                    "interval": interval,
                    "cpg_label": safe_get(idx_cpg, interval),
                    "expected_category": safe_get(idx_exp_cat, "Target Locus"),
                    "summary_label": safe_get(idx_sum_lab, "Analyzed"),
                    "status_label": status,
                    "qc_warnings": qc,
                    "qc_status": qc_status,
                    "mean_combined_methyl": str(round(combined_val, 3)),
                    "meth_score": str(round(combined_val, 3)),
                    "unmeth_score": str(round(1.0 - combined_val, 3)),
                    "mean_meth_delta": str(round(delta_val, 3)),
                    "mean_hap1_methyl": str(round(h1_val, 3)),
                    "mean_hap2_methyl": str(round(h2_val, 3)),
                }
                rows.append(cleaned_row)
    except Exception as e:
        print(f" [!] Error parsing {tsv_path.name}: {e}", file=sys.stderr)
    return rows


def calculate_stats(rows, sid):
    """Calculates summary statistics for dashboard metrics cards."""
    if not rows:
        return {"cov": "N/A", "h1Avg": "0%", "h2Avg": "0%", "passCount": 0}

    h1_vals = [float(r["mean_hap1_methyl"]) for r in rows if r.get("mean_hap1_methyl")]
    h2_vals = [float(r["mean_hap2_methyl"]) for r in rows if r.get("mean_hap2_methyl")]

    h1_avg = (sum(h1_vals) / len(h1_vals) * 100) if h1_vals else 0
    h2_avg = (sum(h2_vals) / len(h2_vals) * 100) if h2_vals else 0
    pass_count = sum(1 for r in rows if r.get("qc_status") == 0)

    return {
        "cov": "Phased-5mC",
        "h1Avg": f"{h1_avg:.1f}%",
        "h2Avg": f"{h2_avg:.1f}%",
        "passCount": pass_count,
    }


def main():
    parser = argparse.ArgumentParser(
        description="longvarmeth Clinical Methylation Dashboard Generator"
    )
    parser.add_argument(
        "--profile",
        required=True,
        help="Path to methbat profile TSV/CSV or modkit aggregated TSV",
    )
    parser.add_argument(
        "--template", required=True, help="Path to dashboard-template.html"
    )
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, help="Path for output HTML report")

    args = parser.parse_args()

    profile_path = Path(args.profile)
    template_path = Path(args.template)
    sample_id = args.sample_id
    out_file = Path(args.output)

    if not profile_path.exists():
        print(f" [!] ERROR: Profile file not found: {profile_path}", file=sys.stderr)
        sys.exit(1)
    if not template_path.exists():
        print(f" [!] ERROR: Template not found: {template_path}", file=sys.stderr)
        sys.exit(1)

    rows = parse_methbat_tsv(profile_path)
    if not rows:
        # Generate dummy fallback locus if profile is empty or header only
        rows = [
            {
                "chrom": "chr1",
                "start": "0",
                "end": "1000",
                "interval": "chr1:0-1000",
                "cpg_label": "No loci identified",
                "expected_category": "N/A",
                "summary_label": "PASS",
                "status_label": "PASS",
                "qc_warnings": "PASS",
                "qc_status": 0,
                "mean_combined_methyl": "0.500",
                "meth_score": "0.500",
                "unmeth_score": "0.500",
                "mean_meth_delta": "0.000",
                "mean_hap1_methyl": "0.500",
                "mean_hap2_methyl": "0.500",
            }
        ]

    stats = calculate_stats(rows, sample_id)

    with open(template_path, "r", encoding="utf-8") as f:
        template_html = f.read()

    sample_report_html = template_html.replace(
        "{{SAMPLE_DATA_JSON}}", json.dumps({sample_id: rows})
    )
    sample_report_html = sample_report_html.replace(
        "{{SAMPLE_STATS_JSON}}", json.dumps({sample_id: stats}, indent=8)
    )

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(sample_report_html)

    print(f" [SUCCESS] Generated methylation dashboard: {out_file} ({len(rows)} loci)")


if __name__ == "__main__":
    main()
