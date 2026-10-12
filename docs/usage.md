# longvarmeth: Usage Guide

This document describes how to configure, parameterize, and run the **longvarmeth** pipeline across **Oxford Nanopore** and **PacBio HiFi** datasets.

---

## Input Samplesheet (`samplesheet.csv`)

The pipeline takes a comma-separated CSV file specifying sample identifiers, data paths, and platform metadata.

### Format

```csv
run_id,input_path,input_type,kit,model,modified_bases_model
sample_ont_raw,/data/raw/pod5/sample_ont_raw,pod5,SQK-NBD114-24,dna_r10.4.1_e8.2_400bps_sup@v5.0.0,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v1
sample_ont_bam,/data/raw/bam/sample_ont.bam,bam,,,
sample_pacbio_hifi,/data/pacbio/m84036_sample.hifi_reads.bam,pacbio_bam,,,
```

### Column Specifications:
* `run_id`: Unique identifier for the sample or barcode.
* `input_path`: Path or URL to a directory of raw `pod5`/`fast5` files, or an existing `.bam` file.
* `input_type`: Input format (`pod5`, `fast5`, `bam`, `pacbio_bam`, `pacbio_hifi`, or `hifi`).
* `kit`: ONT sequencing / barcoding kit code (for host Dorado demultiplexing; optional for PacBio).
* `model`: Dorado basecalling model name or path (optional for PacBio).
* `modified_bases_model`: Dorado modified bases model name (e.g. `5mCG_5hmCG`; optional for PacBio).

---

## Modular Tool Selection

You can easily select your preferred caller stack via CLI flags or in `conf/params.config`:

```bash
--platform        # 'ont' (default), 'pacbio', or 'auto'
--snv_caller      # 'bcftools' (default), 'deepvariant', or 'clair3'
--sv_caller       # 'sniffles' (default), 'severus', or 'pbsv'
--phasing_tool    # 'whatshap' (default for ONT), 'hiphase' (default for PacBio)
--cnv_caller      # 'cnvkit' (default), 'hificnv' (PacBio)
--enable_trgt     # true / false (Tandem Repeat Genotyping via TRGT)
--repeat_bed      # Path to tandem repeat catalog BED file
```

---

## Clinical Governance & PHI Sanitization (from `JD2112/milou`)

To ensure deterministic clinical workflows, audit trails, and strict privacy protections:

```bash
--enable_phi_sanitization true  # Scrub patient IDs and hospital barcodes from BAM/VCF headers
--enable_audit_trail true       # Compute cryptographic SHA-256 manifests for all emitted files
```

When enabled, `<outdir>/audit/` will contain:
* `<sample>_checksums.sha256`: Cryptographic hashes for all published artifacts.
* `<sample>_manifest.json`: Full machine-readable execution record (parameter hashes, git commit, tool versions).

---

## Species Profiles

Preconfigured profiles optimize parameters for specific biology:

* `-profile human`: Diploid (`--ploidy 2`), GRCh38 defaults, WhatsHap/HiPhase phasing enabled.
* `-profile eukaryote`: Diploid ploidy, reference-guided or Hifiasm/Flye assembly with BUSCO lineage evaluation.
* `-profile prokaryote`: Haploid (`--ploidy 1`), single circular reference or Flye/Hifiasm circularization + bacterial 4mC/6mA/5mC motif detection.
* `-profile virus`: Haploid (`--ploidy 1`), quasispecies intra-host diversity calling (relaxed 5% variant frequency cutoffs).

---

## Execution Profiles & Testing

Run pipeline test suites out of the box using public `nf-core/test-datasets`:

```bash
# Nanopore test run:
nextflow run main.nf -profile test_ont,singularity

# PacBio HiFi test run:
nextflow run main.nf -profile test_pacbio,singularity

# Syntax and DAG verification:
nextflow run main.nf -profile test_pacbio -stub
```
