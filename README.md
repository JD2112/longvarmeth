# longvarmeth

[![quindecagon](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/quindecagon.json)](https://github.com/JD2112/quindecagon)
[![nf-core lint](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/nfcore_lint.json)](https://nf-co.re)
[![SBOM](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/syft.json)](https://github.com/anchore/syft)
[![secrets](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/secrets.json)](https://github.com/gitleaks/gitleaks)
[![code style](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/black.json)](https://github.com/psf/black)
[![flake8](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/flake8.json)](https://flake8.pycqa.org)
[![R audit](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/r_audit.json)](https://github.com/r-lib/lintr)

A modular, clinical-grade Nextflow DSL2 platform for single-molecule long-read genomics and epigenomics (**Oxford Nanopore & PacBio HiFi**). It supports host GPU-accelerated basecalling/demux, comprehensive QC, *de novo* assembly & polishing, reference alignment, epigenetic tag preservation (`MM`/`ML`), multi-tool variant calling (SNV/InDel/SV/CNV/STR), haplotype co-phasing, strict clinical governance & PHI sanitization, and interactive visual reporting (JBrowse 2 / Methylartist / MultiQC).

---

## Metro Architecture

![longvarmeth Metro Diagram](./docs/pipeline.png)

---

## Key Features

* **Multi-Platform Long-Read Ingestion**:
  * **Oxford Nanopore (ONT)**: Host GPU-accelerated basecalling & demultiplexing via ONT **Dorado** from raw `pod5`/`fast5` or pre-basecalled BAMs.
  * **PacBio HiFi / Revio**: Native PacBio BAM ingestion preserving on-instrument kinetic modification tags (`MM`/`ML`).
* **Clinical Governance & Forensic Audit Trail (Inspired by `JD2112/milou`)**:
  * **PHI Sanitization**: Automatically scrubs patient identifiers, hospital run barcodes, and sensitive comments from BAM and VCF headers before publication.
  * **Cryptographic Provenance**: Emits SHA-256 manifests (`<sample>_checksums.sha256`) and structured execution audit logs (`<sample>_manifest.json`).
* **Quality Control & Assembly**:
  * Multi-tool read assessment with **FastQC**, **Cramino**, **Mosdepth**, and **NanoPlot**.
  * Dual *de novo* assembly (**Flye** for ONT; **Hifiasm** for PacBio), consensus polishing (**Medaka** / **Racon**), and lineage completeness evaluation (**QUAST**, **BUSCO**).
* **Reference Alignment**:
  * ONT alignment via **Minimap2** (`-ax map-ont -y`).
  * PacBio HiFi alignment via **pbmm2** (`--preset HIFI`) or Minimap2 (`-ax map-hifi -y`), preserving methylation kinetics.
* **Variant Calling & Phasing Matrix**:
  * **SNVs / InDels**: High-accuracy calling via **DeepVariant**, **Clair3**, or **BCFtools**.
  * **Structural Variants (SVs)**: High-sensitivity SV detection with **Sniffles2** and **Severus**.
  * **Tandem Repeat Expansions (STRs)**: Pathogenic repeat genotyping using **TRGT**.
  * **Copy Number Alterations (CNVs)**: Whole-genome copy number profiling via **HiFiCNV** (PacBio) and **CNVkit** (ONT / cross-platform).
  * **Haplotype Phasing**: Read haplotagging and co-phasing of SNVs, SVs, and STRs via **HiPhase** (PacBio) or **WhatsHap** (ONT).
* **Epigenomics & Base Modifications**:
  * Whole-genome 5mC / 5hmC pileups with ONT **Modkit**.
  * PacBio-native CpG aggregation via **pbcpgtools**.
  * Differentially methylated region (DMR) calling with **methbat** and `modkit dmr`.
  * Bacterial restriction-modification motif identification (4mC, 6mA, 5mC).
* **Interactive Visualization**:
  * Automatically generates **JBrowse 2** interactive browser tracks (with colored modification tags and phased tracks).
  * Automated target locus plotting with **Methylartist**.
* **Four Dedicated Organism Profiles**:
  * `-profile human`: Diploid clinical/germline WGS/WES with variant annotation and phasing.
  * `-profile eukaryote`: Diploid non-model animals, plants, and fungi with assembly and BUSCO QC.
  * `-profile prokaryote`: Haploid bacteria/archaea single-contig circularization and epigenetic motif extraction.
  * `-profile virus`: Viral quasispecies and intra-host variant calling (SARS-CoV-2, HIV, phages, amplicons).

---

## Quickstart

Run the pipeline using a specific platform and species profile:

### 1. Oxford Nanopore Reference-Based Run (Human / Eukaryote)
```bash
nextflow run main.nf \
  -profile human,singularity \
  --input samplesheet.csv \
  --ref_nuc /path/to/GRCh38.fa \
  --outdir results_ont_human
```

### 2. PacBio HiFi Reference-Based Run
```bash
nextflow run main.nf \
  -profile human,singularity \
  --platform pacbio \
  --input samplesheet_pacbio.csv \
  --ref_nuc /path/to/GRCh38.fa \
  --enable_trgt \
  --repeat_bed /path/to/pathogenic_repeats.bed \
  --outdir results_pacbio_human
```

### 3. Bacterial / Prokaryotic Assembly & Methylation
```bash
nextflow run main.nf \
  -profile prokaryote,singularity \
  --input samplesheet.csv \
  --run_type de_novo_assembly \
  --outdir results_bacteria
```

### 4. Viral Genome Analysis
```bash
nextflow run main.nf \
  -profile virus,singularity \
  --input samplesheet_viral.csv \
  --ref_nuc /path/to/viral_ref.fasta \
  --outdir results_viral
```

---

## Testing

The pipeline includes built-in test profiles utilizing public `nf-core/test-datasets`:

```bash
# Test Nanopore workflow
nextflow run main.nf -profile test_ont,singularity

# Test PacBio HiFi workflow
nextflow run main.nf -profile test_pacbio,singularity

# Fast topology/syntax check without downloading containers
nextflow run main.nf -profile test_pacbio -stub
```

---

## Species Profiles

| Profile | Target Organism Type | Key Settings |
| :--- | :--- | :--- |
| `-profile human` | *Homo sapiens* | Diploid (`--ploidy 2`), hg38 defaults, WhatsHap/HiPhase phasing, VEP/SnpEff enabled |
| `-profile eukaryote` | Non-Model Eukaryotes | Diploid ploidy, reference-guided or Hifiasm/Flye assembly with BUSCO lineage QC |
| `-profile prokaryote`| Bacteria & Archaea | Haploid (`--ploidy 1`), circularized assembly, bacterial 4mC/6mA/5mC motif detection |
| `-profile virus` | Viruses (Amplicons / Quasispecies) | Haploid (`--ploidy 1`), relaxed 5% frequency cutoffs for intra-host diversity |

---

## Output Structure

Detailed outputs are published in `--outdir` (default: `results/`):

```
results/
├── qc/                       # FastQC, Cramino, NanoPlot, Mosdepth reports
├── assembly/                 # Flye / Hifiasm assemblies, Medaka polished contigs, QUAST & BUSCO
├── alignments/               # Coordinate-sorted BAMs, MM/ML tagged BAMs, haplotagged BAMs
├── variants/                 # SNV VCFs, Sniffles2 SVs, TRGT STRs, HiFiCNV / CNVkit CNVs
├── methylation/              # Modkit 5mC/5hmC pileup BED files, pbcpgtools scores, DMR tables
├── audit/                    # Cryptographic SHA-256 manifests and forensic provenance records
├── jbrowse/                  # Interactive JBrowse 2 web workspace directory
└── multiqc/                  # Consolidated multi-sample MultiQC report
```
See [docs/output.md](docs/output.md) for full descriptions.
