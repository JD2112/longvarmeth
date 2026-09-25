# longvarmeth: Output Guide

This document describes the directory structure and files generated in the `--outdir` (default: `results/`).

```
results/
├── qc/
│   ├── fastqc/               # Raw read FastQC HTML reports and ZIP archives
│   ├── nanoplot/             # Read length, yield, and quality distribution plots
│   ├── cramino/              # Alignment and yield QC metrics for BAM files
│   ├── filtlong/             # Size and quality filtered FASTQ files (assembly mode)
│   ├── quast/                # Assembly contiguity and statistics reports
│   └── busco/                # Assembly gene completeness benchmark reports
├── assembly/
│   ├── flye/                 # Flye de novo assembly FASTA, GFA graphs, and logs
│   └── medaka/               # Medaka consensus polished FASTA files
├── alignments/
│   ├── sorted.bam (.bai)     # Coordinate-sorted genomic alignments
│   ├── tagged/               # BAM files with restored methylation MM/ML tags
│   └── haplotagged/          # WhatsHap phased and haplotagged BAM files (HP:i:1, HP:i:2)
├── variants/
│   ├── raw/                  # Per-chromosome raw BCF/VCF files
│   ├── <sample>.vcf.gz       # Merged, genome-wide small variant calls (SNVs/InDels)
│   ├── phased/               # WhatsHap phased VCF files
│   ├── svs/                  # Sniffles2 structural variant calls (.vcf.gz)
│   └── cnvs/                 # CNVkit copy number segments (.cns) and diagram plots (.pdf)
├── methylation/
│   └── <sample>_mods.bed     # Modkit 5mC/5hmC methylation frequency pileup BED files
├── jbrowse/
│   └── <sample>/             # Standalone interactive JBrowse 2 web instances
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

## Detailed Directory Descriptions

### `qc/`
* **`fastqc/`**: Standard FastQC results per barcode/sample.
* **`nanoplot/`**: Yield vs. length scatter plots, read length histograms, and logarithmic quality distributions.
* **`quast/`**: Assembly metrics including N50, L50, total length, and GC%.
* **`busco/`**: Assessment of single-copy orthologs present in polished assemblies.

### `alignments/`
* **`alignments/tagged/`**: Contains `<sample>.mmtags.sorted.bam`. These BAM files carry the original `MM` and `ML` SAM tags representing 5mC or 6mA modifications preserved through Minimap2 alignment.
* **`alignments/haplotagged/`**: Generated when phasing is enabled. Reads are tagged with `HP:i:1` or `HP:i:2` for allele-specific visualization in genome browsers.

### `variants/`
* **`<sample>.vcf.gz`**: High-confidence small variants (SNVs and InDels) filtered by quality and coverage.
* **`svs/<sample>_svs.vcf.gz`**: Structural variants (deletions, duplications, inversions, translocations) called by Sniffles2.
* **`cnvs/`**: Copy number alteration profiles and whole-genome CNV diagrams produced by CNVkit.

### `methylation/`
* **`<sample>_mods.bed`**: BedMethyl format reporting total coverage, modified base counts, and methylation percentages at every CpG site across the genome.

### `jbrowse/`
* Self-contained JBrowse 2 instances with the reference FASTA, tagged/haplotagged BAMs, and VCF variant tracks pre-configured. Can be served via `python3 -m http.server 8000`.

### `multiqc/`
* **`multiqc_report.html`**: Complete quality summary combining FastQC, NanoPlot, Cramino, QUAST, BUSCO, and software versions.
