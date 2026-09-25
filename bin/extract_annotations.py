#!/usr/bin/env python3


# Here is a detailed breakdown of the column headers and what the values in your output mean for researchers:

# ### 1. Columns and Values Explained (from your example)

# | Column Header | Value in Example | Explanation |
# | :--- | :--- | :--- |
# | **Sample** | `sample_bc05.mmtags...` | The name of the sample (Barcode 05). |
# | **Gene** | `MEFV` | The target gene. |
# | **Variant (Protein)** | *(Blank)* | The amino acid change. It is **blank** here because this variant is non-coding (located downstream of the coding sequence), so there is no protein change. |
# | **Variant (cDNA)** | `c.*4805G>A` | The cDNA nucleotide position change. The `*` prefix means this is in the 3' untranslated (UTR) or downstream region. |
# | **Genotype** | `Homozygous Alt` | Both chromosome copies have the alternative allele (mutant genotype). |
# | **Site_Methylation (%)** | `N/A` | The CpG methylation level at this exact site. It is `N/A` because this variant position is not a CpG site. |
# | **Site_Methylation_Coverage** | `N/A` | The read depth at this exact site for methylation basecalling. |
# | **Effect** | `downstream_gene_variant` | The SnpEff predicted sequence effect (e.g. missense, synonymous, downstream). |
# | **Impact** | `MODIFIER` | SnpEff's predicted functional severity (HIGH, MODERATE, LOW, or MODIFIER). |
# | **Chromosome** | `chr16` | The chromosome location. |
# | **Position** | `3238336` | The exact 1-based genomic coordinate (hg38). |
# | **Ref** | `C` | The reference nucleotide allele. |
# | **Alt** | `T` | The mutant/alternative nucleotide allele. |
# | **rsID** | `.` | The dbSNP rsID database identifier. `.` means this is a novel or undocumented variant. |
# | **GQ** | `127` | Genotype Quality. A score of `127` indicates an extremely high-confidence variant call. |

# ---

# ### 2. Header and Column Details for both output files

# #### File A: `target_gene_variants.csv`
# This file lists every sequence variant detected within the boundaries of the five target genes:
# 1. **Sample**: The identifier of the sequencing sample.
# 2. **Gene**: Gene name (e.g. `MEFV`).
# 3. **Variant (Protein)**: The HGVS protein notation (e.g. `p.M694V`).
# 4. **Variant (cDNA)**: The HGVS cDNA notation (e.g. `c.2080A>G`).
# 5. **Genotype**: The genotype of the sample (e.g., `Heterozygous`, `Homozygous Alt`).
# 6. **Site_Methylation (%)**: Percentage of reads methylated at this exact single base coordinate.
# 7. **Site_Methylation_Coverage**: Number of reads covering this site for methylation assessment.
# 8. **Effect**: Sequence ontology consequence (e.g. `missense_variant`, `synonymous_variant`, `intron_variant`).
# 9. **Impact**: Classification of effect severity.
# 10. **Chromosome / Position / Ref / Alt / rsID**: Genomic coordinates and alleles.
# 11. **GQ**: Phred-scaled confidence score that the called genotype is correct.

# #### File B: `target_gene_methylation_summary.csv`
# This file summarizes the epigenetic methylation level across the entire region of each gene:
# 1. **Sample**: The identifier of the sequencing sample.
# 2. **Gene**: Target gene.
# 3. **Chrom / Start / End**: The hg38 coordinate boundaries used to define the gene.
# 4. **CpG_Sites_Detected**: The total number of valid CpG loci covered by reads within this gene.
# 5. **Average_Methylation (%)**: The mean methylation percentage across all detected CpG sites in this gene. *(Perfect for comparing promoter hyper- or hypo-methylation across patients)*.




import gzip
import glob
import re
import os
import sys

# Target genes and their hg38 coordinates
GENE_REGIONS = {
    "MEFV": {"chrom": "chr16", "start": 3240392, "end": 3256954},
    "NLRP3": {"chrom": "chr1", "start": 247332331, "end": 247449668},
    "TNFRSF1A": {"chrom": "chr12", "start": 6328757, "end": 6342114},
    "MVK": {"chrom": "chr12", "start": 109573255, "end": 109598125},
    "TNFAIP3": {"chrom": "chr6", "start": 137866383, "end": 137885836}
}

TARGET_GENES = set(GENE_REGIONS.keys())

def load_methylation_data(bed_path):
    """
    Parses a bedMethyl file and returns:
    1. A dictionary of exact site methylation: (chrom, pos_1based) -> (percent, coverage)
    2. A list of all CpG site records: (chrom, pos_1based, percent, coverage)
    """
    site_dict = {}
    all_sites = []
    
    if not os.path.exists(bed_path):
        return site_dict, all_sites
        
    open_func = gzip.open if bed_path.endswith('.gz') else open
    mode = 'rt' if bed_path.endswith('.gz') else 'r'
    
    try:
        with open_func(bed_path, mode) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                chrom = fields[0]
                pos_1based = int(fields[2])  # 1-based end coordinate
                coverage = int(fields[9])
                percentage = float(fields[10])
                
                site_dict[(chrom, pos_1based)] = (percentage, coverage)
                all_sites.append((chrom, pos_1based, percentage, coverage))
    except Exception as e:
        print(f"Warning: Failed to parse methylation file {bed_path}: {e}", file=sys.stderr)
        
    return site_dict, all_sites

def calculate_gene_methylation(all_sites, sample_id):
    """
    Computes average methylation stats for the target gene regions.
    """
    stats = []
    for gene, region in GENE_REGIONS.items():
        chrom = region["chrom"]
        start = region["start"]
        end = region["end"]
        
        # Filter CpG sites in the gene region
        gene_sites = [
            percent for c, pos, percent, cov in all_sites 
            if c == chrom and start <= pos <= end
        ]
        
        if gene_sites:
            avg_meth = sum(gene_sites) / len(gene_sites)
            num_sites = len(gene_sites)
        else:
            avg_meth = 0.0
            num_sites = 0
            
        stats.append({
            "Sample": sample_id,
            "Gene": gene,
            "Chrom": chrom,
            "Start": start,
            "End": end,
            "CpG_Sites_Detected": num_sites,
            "Average_Methylation (%)": f"{avg_meth:.2f}" if num_sites > 0 else "N/A"
        })
    return stats

def parse_vcf(vcf_path, meth_site_dict):
    sample_id = os.path.basename(vcf_path).split('.')[0]
    results = []
    
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
    with open_func(vcf_path, mode) as f:
        samples_headers = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                samples_headers = headers[9:]
                continue
            
            # Parse variant lines
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            rsid = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            format_field = fields[8].split(':')
            samples_data = fields[9:]
            
            # Find genotype and quality indices
            gt_idx = format_field.index('GT') if 'GT' in format_field else -1
            gq_idx = format_field.index('GQ') if 'GQ' in format_field else -1
            
            # Extract SnpEff annotations (ANN field)
            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue
            
            annotations = ann_match.group(1).split(',')
            for ann in annotations:
                ann_fields = ann.split('|')
                if len(ann_fields) < 11:
                    continue
                
                gene = ann_fields[3]
                if gene in TARGET_GENES:
                    effect = ann_fields[1]
                    impact = ann_fields[2]
                    hgvs_c = ann_fields[9]
                    hgvs_p = ann_fields[10]
                    
                    # Look up methylation info at this exact position
                    meth_percent = "N/A"
                    meth_cov = "N/A"
                    if (chrom, pos) in meth_site_dict:
                        meth_percent, meth_cov = meth_site_dict[(chrom, pos)]
                        meth_percent = f"{meth_percent:.2f}%"
                    
                    # Extract sample genotypes
                    for idx, sample_val in enumerate(samples_data):
                        sample_name = samples_headers[idx] if idx < len(samples_headers) else sample_id
                        vals = sample_val.split(':')
                        gt = vals[gt_idx] if gt_idx < len(vals) else './.'
                        gq = vals[gq_idx] if gq_idx < len(vals) else '.'
                        
                        gt_desc = "Unknown"
                        if gt in ["0/1", "0|1", "1|0"]:
                            gt_desc = "Heterozygous"
                        elif gt in ["1/1", "1|1"]:
                            gt_desc = "Homozygous Alt"
                        elif gt in ["0/0", "0|0"]:
                            continue  # Skip reference homozygotes
                            
                        results.append({
                            "Sample": sample_name,
                            "Gene": gene,
                            "Variant (Protein)": hgvs_p,
                            "Variant (cDNA)": hgvs_c,
                            "Genotype": gt_desc,
                            "Site_Methylation (%)": meth_percent,
                            "Site_Methylation_Coverage": meth_cov,
                            "Effect": effect,
                            "Impact": impact,
                            "Chromosome": chrom,
                            "Position": pos,
                            "Ref": ref,
                            "Alt": alt,
                            "rsID": rsid,
                            "GQ": gq
                        })
                    break  # Output one annotation per variant site
                    
    return results

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract target gene variants and methylation summaries.")
    parser.add_argument("--input-dir", default="results", help="Nextflow pipeline output directory (defaults to 'results')")
    args = parser.parse_args()

    vcf_dir = os.path.join(args.input_dir, "annotated_variants")
    if not os.path.exists(vcf_dir):
        # Fallback to checking the variants/ directory directly in newer versions of the pipeline
        vcf_dir = os.path.join(args.input_dir, "variants")

    vcf_pattern = os.path.join(vcf_dir, "*.snpeff.vcf.gz")
    vcf_files = glob.glob(vcf_pattern)
    
    if not vcf_files:
        # Fallback to checking for standard annotation or unannotated files
        vcf_pattern = os.path.join(vcf_dir, "*.vcf.gz")
        vcf_files = glob.glob(vcf_pattern)
        
    if not vcf_files:
        print(f"No VCF files found in {vcf_dir} matching pattern *.vcf.gz", file=sys.stderr)
        sys.exit(1)
        
    all_variants = []
    all_gene_methylation = []
    
    for vcf_file in vcf_files:
        filename = os.path.basename(vcf_file)
        sample_id = filename.split('.')[0].replace('.snpeff', '')
        
        # Locate corresponding methylation bed file
        # Check both .bed and .bed.gz
        meth_bed = os.path.join(args.input_dir, "methylation", f"{sample_id}_mods.bed")
        if not os.path.exists(meth_bed):
            meth_bed_gz = f"{meth_bed}.gz"
            if os.path.exists(meth_bed_gz):
                meth_bed = meth_bed_gz
                
        # Load methylation data
        meth_site_dict, all_meth_sites = load_methylation_data(meth_bed)
        
        # 1. Extract variants and annotate with site methylation
        all_variants.extend(parse_vcf(vcf_file, meth_site_dict))
        
        # 2. Extract gene-level methylation stats
        all_gene_methylation.extend(calculate_gene_methylation(all_meth_sites, sample_id))
        
    # Determine the output prefix from the input directory name
    dir_name = os.path.basename(os.path.normpath(args.input_dir))
    
    # Write variant table
    variant_outfile = f"{dir_name}_target_gene_variants.csv"
    var_columns = [
        "Sample", "Gene", "Variant (Protein)", "Variant (cDNA)", "Genotype", 
        "Site_Methylation (%)", "Site_Methylation_Coverage", "Effect", "Impact", 
        "Chromosome", "Position", "Ref", "Alt", "rsID", "GQ"
    ]
    with open(variant_outfile, "w") as f:
        f.write(",".join(var_columns) + "\n")
        all_variants.sort(key=lambda x: (x["Sample"], x["Gene"], x["Chromosome"], int(x["Position"])))
        for r in all_variants:
            row = [str(r[col]) for col in var_columns]
            f.write(",".join(row) + "\n")
    print(f"Successfully wrote variant stats to: {variant_outfile}")
            
    # Write gene-level methylation summary table
    meth_outfile = f"{dir_name}_target_gene_methylation_summary.csv"
    meth_columns = ["Sample", "Gene", "Chrom", "Start", "End", "CpG_Sites_Detected", "Average_Methylation (%)"]
    with open(meth_outfile, "w") as f:
        f.write(",".join(meth_columns) + "\n")
        all_gene_methylation.sort(key=lambda x: (x["Sample"], x["Gene"]))
        for r in all_gene_methylation:
            row = [str(r[col]) for col in meth_columns]
            f.write(",".join(row) + "\n")
    print(f"Successfully wrote gene methylation stats to: {meth_outfile}")

if __name__ == "__main__":
    main()
