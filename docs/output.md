# longvarmeth: Output Guide

This document describes the directory structure and files generated in `--outdir` (default: `results/`).

```
results/
├── qc/
│   ├── fastqc/               # Raw read FastQC HTML reports and ZIP archives
│   ├── nanoplot/             # Read length, yield, and quality distribution plots
│   ├── cramino/              # Alignment and yield QC metrics for BAM files
│   ├── mosdepth/             # High-speed whole-genome coverage and depth statistics
│   ├── filtlong/             # Size and quality filtered FASTQ files (assembly mode)
│   ├── quast/                # Assembly contiguity and statistics reports
│   └── busco/                # Assembly gene completeness benchmark reports
├── assembly/
│   ├── flye/                 # Flye de novo assembly FASTA, GFA graphs, and logs
│   ├── hifiasm/              # Hifiasm phased diploid assemblies (hap1, hap2, primary)
│   └── medaka/               # Medaka consensus polished FASTA files
├── alignments/
│   ├── pbmm2/                # Coordinate-sorted PacBio HiFi alignments (.bam, .bai)
│   ├── sorted.bam (.bai)     # Coordinate-sorted genomic alignments (Minimap2)
│   ├── tagged/               # BAM files with restored methylation MM/ML tags
│   └── haplotagged/          # Phased and haplotagged BAM files (HP:i:1, HP:i:2)
├── variants/
│   ├── raw/                  # Per-chromosome raw BCF/VCF files
│   ├── <sample>.vcf.gz       # Small variant calls (SNVs/InDels) via BCFtools / DeepVariant / Clair3
│   ├── phased/               # Phased VCF files (WhatsHap / HiPhase)
│   ├── svs/                  # Sniffles2 / Severus structural variant calls (.vcf.gz)
│   ├── trgt/                 # TRGT tandem repeat genotypes and spanning read BAMs
│   └── cnvs/                 # HiFiCNV / CNVkit copy number calls and diagrams
├── methylation/
│   ├── <sample>_mods.bed     # Modkit 5mC/5hmC methylation frequency pileup BED files
│   ├── pbcpgtools/           # PacBio CpG aggregated scores
│   └── dmr/                  # Phased differentially methylated regions (DMRs)
├── audit/
│   ├── <sample>_checksums.sha256 # Cryptographic SHA-256 manifest of all published artifacts
│   └── <sample>_manifest.json    # Machine-readable execution provenance log (tools, seeds, parameters)
├── jbrowse/
│   └── <sample>/             # Standalone interactive JBrowse 2 web instances
├── reports/
│   └── methylation/          # Interactive clinical methylation dashboards (AG-Grid, IGV-Web, gnomAD links)
├── multiqc/
│   ├── multiqc_report.html   # Consolidated interactive QC summary report
│   └── multiqc_data/         # Aggregated metrics tables and logs
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    ├── pipeline_dag.html
    └── software_versions_mqc.yaml
```

---

## Key Directories Explained

### `audit/` (Clinical Governance)
Generated when `--enable_audit_trail true` or `--enable_phi_sanitization true` is set:
* **`<sample>_checksums.sha256`**: Exact cryptographic digests of every final BAM, VCF, and BED file for data integrity verification.
* **`<sample>_manifest.json`**: Forensic metadata log documenting the exact Nextflow commit, container hashes, and parameter states.

### `variants/trgt/`
Emitted when `--enable_trgt true` is specified:
* **`<sample>.trgt.vcf.gz`**: Genotyped repeat alleles with repeat motif counts, methylation levels per allele, and confidence intervals.
* **`<sample>.spanning.bam`**: Reads spanning repeat loci tagged with repeat-allele assignments for visual inspection.

### `alignments/`
* **`*.bam` / `*.bai`**: Coordinate-sorted, indexed alignments preserving raw sequencing kinetic tags (`MM`, `ML`).
* **`haplotagged/*.bam`**: Read alignments with added haplotype assignment tags (`HP:i:1`, `HP:i:2`) for allele-specific visualization in IGV and JBrowse 2.

### `reports/methylation/` (Clinical Methylation Dashboard)
Generated when `--enable_methylation_dashboard true` is set (enabled automatically with `-profile human`):
* **`<sample>_methylation_dashboard.html`**: Self-contained clinical web interface adapted from `gms-nallo-methylation-dashboard`. Features an interactive AG-Grid triage table (DMRs, expected categories, QC filters), embedded lazy-loading IGV-Web showing haplotagged phased reads and base modifications (5mC/5hmC), and instant one-click hyperlinks to gnomAD v4 and NCBI GDV.

* **`alignments/tagged/`**: Contains `<sample>.sanitized.bam` or `<sample>.mmtags.sorted.bam`. These BAM files preserve the original `MM` and `ML` SAM tags representing 5mC, 5hmC, or 6mA modifications without stripping downstream headers.
