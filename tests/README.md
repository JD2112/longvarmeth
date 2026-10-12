# Pipeline Validation & Testing Guide

This directory contains test configurations, sample sheets, and instructions for running end-to-end and stub tests across Oxford Nanopore and PacBio platforms on your server (`fraka`).

---

## 1. Nanopore Test

### Stub Dry-Run (Fast Syntax & Topology Check):
```bash
nextflow run main.nf -profile test_ont -stub
# or:
nextflow run main.nf -profile test -stub
```

### Full Run on Server:
```bash
nextflow run main.nf -profile test_ont,fraka
```

---

## 2. PacBio HiFi Test (using nf-core test datasets)

Tests the pipeline with official PacBio HiFi test data (`m84036_240201_205955_s2.hifi_reads.bc2001_subsample.bam` + `C9ORF72-12.fa`) hosted on `nf-core/test-datasets`:

### Stub Dry-Run:
```bash
nextflow run main.nf -profile test_pacbio -stub
```

### Full Run with Containerized Execution on Server:
```bash
nextflow run main.nf -profile test_pacbio,singularity
# or with fraka config:
nextflow run main.nf -profile test_pacbio,fraka
```

---

## 3. Viral Genome Test (Profile `virus`)

```bash
nextflow run main.nf \
    -profile virus,fraka \
    --input samplesheet_viral.csv \
    --ref_nuc https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta \
    --outdir results_viral
```

---

## 4. Clinical Governance & Audit Outputs

When `--enable_phi_sanitization` and `--enable_audit_trail` are active:
- Sanitized BAMs are emitted without clinical identifying comments/headers.
- Cryptographic SHA-256 manifests (`<sample>_checksums.sha256`) and provenance manifests (`<sample>_manifest.json`) are automatically written to `<outdir>/audit/`.
