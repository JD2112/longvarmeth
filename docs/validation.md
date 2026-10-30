# Analytical Validation & Truth-Set Benchmarking Guide

## 1. Intended Use and Validation Scope

**LongVarMeth is intended for Research Use Only (RUO) and clinical evaluation.**  
It provides a rigorous bioinformatics pipeline incorporating clinical governance controls (cryptographic SHA-256 manifests, metadata and PHI scrubbing, deterministic subject de-identification). However, before using pipeline results for medical diagnosis, clinical management, or patient reporting, clinical laboratories must independently establish analytical validity (accuracy, precision, analytical sensitivity, analytical specificity, and reportable range) in compliance with applicable standards (e.g. CLIA/CAP, ISO 15189, or regional IVD directives).

---

## 2. GIAB HG002 Truth-Set Benchmarking Architecture

To provide transparent, empirical verification of variant calling accuracy, LongVarMeth includes an automated benchmarking subworkflow based on:
1. **`hap.py` (v0.3.15)**: The industry-standard variant comparison tool developed by Illumina and the Global Alliance for Genomics and Health (GA4GH) Benchmarking Team for small variants (SNVs and Indels).
2. **`Truvari` (v4.2.2)**: The standard benchmarking tool for structural variants ($\ge 50$ bp) endorsed by the Genome in a Bottle (GIAB) Consortium.

### Benchmark Data Reference
- **Sample**: GIAB HG002 / NA24385 (Ashkenazim Son)
- **Small Variants Truth Set**: NIST GIAB v4.2.1 Small Variant Benchmark
- **Structural Variants Truth Set**: NIST GIAB v0.6 / v1.0 Tier 1 Structural Variant Benchmark
- **High-Confidence Regions**: Stratified BED files defining high-confidence genomic intervals.

---

## 3. Running the Benchmarking Subworkflow

### Using the Preconfigured GIAB Profile (`-profile test_giab`)
```bash
nextflow run main.nf \
  -profile test_giab \
  -stub
```

### Running with Real GIAB HG002 Data
To evaluate an actual GIAB HG002 run against the NIST truth set:

```bash
nextflow run main.nf \
  -profile human,singularity,version_lock \
  --platform ont \
  --input samplesheet_hg002.csv \
  --ref_nuc /path/to/GRCh38_no_alt.fa \
  --run_benchmarking true \
  --benchmark_truth_vcf /data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --benchmark_truth_bed /data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --benchmark_sv_truth_vcf /data/truth/HG002_SVs_Tier1_v0.6.vcf.gz \
  --benchmark_sv_truth_bed /data/truth/HG002_SVs_Tier1_v0.6.bed \
  --outdir results_hg002_benchmark
```

---

## 4. Stratified Performance Expectations

The table below summarizes expected baseline concordance metrics across supported long-read callers on GIAB HG002:

| Variant Type | Tool / Caller | Sequencing Platform | Coverage | High-Confidence Regions (F1) | Difficult / Homopolymer Regions (F1) |
| :--- | :--- | :--- | :--- | :---: | :---: |
| **SNV** | `DeepVariant v1.6.1` (`WGS_PACBIO`) | PacBio HiFi (Revio) | 30x | **> 0.999** | **> 0.995** |
| **SNV** | `DeepVariant v1.6.1` (`ONT_R104`) | ONT R10.4.1 (Simplex) | 35x | **> 0.996** | **> 0.985** |
| **SNV** | `Clair3 v1.0.10` | ONT R10.4.1 (Simplex) | 35x | **> 0.995** | **> 0.980** |
| **Indel** | `DeepVariant v1.6.1` (`WGS_PACBIO`) | PacBio HiFi (Revio) | 30x | **> 0.994** | **> 0.982** |
| **Indel** | `DeepVariant v1.6.1` (`ONT_R104`) | ONT R10.4.1 (Simplex) | 35x | **> 0.880** | **0.750 – 0.840** |
| **Indel** | `Clair3 v1.0.10` | ONT R10.4.1 (Simplex) | 35x | **> 0.875** | **0.740 – 0.830** |
| **SV ($\ge 50$ bp)** | `Sniffles2 v2.8.0` | PacBio HiFi / ONT | $\ge 25$x | **> 0.945** | **0.880 – 0.920** |

> [!WARNING]
> **Homopolymer Indel Limitation on Nanopore Data**:  
> In simplex Oxford Nanopore R10.4.1 sequencing, indel recall in homopolymer tracts ($N \ge 6$) is lower than in non-homopolymer regions due to unresolvable ionic current saturation. High-consequence clinical indel findings in homopolymers must be orthogonally confirmed (e.g. Sanger sequencing or PacBio HiFi).

---

## 5. MultiQC Execution Status Matrix & Benchmark Reporting

Upon completion of a benchmark-enabled run, the MultiQC report will feature:
1. **Pipeline Execution Status Matrix**: Details the exact status of all 12 stages (`COMPLETED`, `SKIPPED`, `FAILED`).
2. **Hap.py Truth-Set Benchmark Table**: Reports SNP and Indel recall, precision, and F1-score.
3. **Truvari Structural Variant Benchmark Table**: Reports SV true positives, false positives, false negatives, recall, precision, and F1-score.
4. **Audit Manifest Summary**: Provides cryptographic SHA-256 hashes for all input files, reference fasta, configuration, and generated benchmark outputs.

