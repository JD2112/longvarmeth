# longvarmeth: Usage Guide

This document describes how to configure, parameterize, and run the **longvarmeth** pipeline.

## Input Samplesheet (`samplesheet.csv`)

The pipeline takes a comma-separated CSV file specifying runs and raw or aligned data.

### Format

```csv
run_id,input_path,input_type,kit,model,modified_bases_model
sample_bc01,/data/raw/pod5/sample_bc01,pod5,SQK-NBD114-24,dna_r10.4.1_e8.2_400bps_sup@v5.0.0,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v1
sample_bc02,/data/raw/bam/sample_bc02.bam,bam,,,
```

### Columns:
* `run_id`: Unique identifier for the sample / barcode.
* `input_path`: Path to a directory of `pod5` / `fast5` files, or an existing `.bam` file.
* `input_type`: Input format (`pod5`, `fast5`, or `bam`).
* `kit`: ONT sequencing / barcoding kit code (for host Dorado demultiplexing).
* `model`: Dorado basecalling model name or path.
* `modified_bases_model`: Dorado modified bases model name (e.g., `5mCG_5hmCG`).

---

## Analysis Run Types

### 1. Reference-Based Mode (`--run_type reference_based`)
Aligns long reads against reference genomes, transfers base modification tags (`MM`/`ML`), calls small variants (SNVs/indels), structural variants (Sniffles), copy number alterations (CNVkit), performs phasing (WhatsHap), and aggregates 5mC methylation pileups (Modkit).

```bash
nextflow run main.nf \
  -profile human,singularity \
  --input samplesheet.csv \
  --run_type reference_based \
  --ref_nuc /path/to/hg38.fa \
  --outdir results_human
```

### 2. De Novo Assembly Mode (`--run_type de_novo_assembly`)
Filters reads using *Filtlong*, generates high-accuracy long-read assemblies with *Flye*, performs consensus polishing with *Medaka*, and runs assembly QC using *QUAST* and *BUSCO*.

```bash
nextflow run main.nf \
  -profile prokaryote,singularity \
  --input samplesheet.csv \
  --run_type de_novo_assembly \
  --outdir results_bacterial_assembly
```

---

## Species Profiles (`-profile <species>`)

Predefined profiles set ploidy, reference genomes, and calling parameters:

* `-profile human`: Diploid (`--ploidy 2`), hg38 reference defaults, WhatsHap phasing enabled.
* `-profile eukaryote`: Diploid (`--ploidy 2`), general eukaryote parameters.
* `-profile prokaryote`: Haploid (`--ploidy 1`), `bacteria_odb10` BUSCO lineage.
* `-profile virus`: Haploid (`--ploidy 1`), small genome parameters.

---

## Execution Engine Profiles

Select your container engine:
* `-profile singularity`: Recommended on shared HPC systems (bind mounts `/data`).
* `-profile apptainer`: Recommended on modern Linux HPC systems.
* `-profile docker`: Recommended for local workstations.
* `-profile conda`: Fallback for local Conda environments.

---

## Key Parameters Summary

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input` | String | `samplesheet.csv` | Input CSV samplesheet |
| `--outdir` | String | `results` | Output directory for all results |
| `--species` | String | `eukaryote` | Species profile (`human`, `eukaryote`, `prokaryote`, `virus`) |
| `--run_type` | String | `reference_based` | `reference_based` or `de_novo_assembly` |
| `--ref_nuc` | Path | `""` | Primary nuclear / chromosomal FASTA reference |
| `--enable_phasing` | Boolean | `false` | Enable WhatsHap variant phasing and haplotagging |
| `--busco_lineage` | String | `eukaryota` | Lineage dataset for BUSCO QC (`bacteria_odb10`, `fungi_odb10`, etc.) |
| `--medaka_model` | String | `r1041_e82_400bps_hac_v4.0.0` | Model for Medaka consensus polishing |
| `--filt_min_len` | Integer | `1000` | Minimum read length for Filtlong filtering (bp) |
| `--filt_keep_pct` | Integer | `90` | Best % bases kept by Filtlong |
