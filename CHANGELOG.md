# Changelog

All notable changes to the **longvarmeth** pipeline are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.2.1] - 2026-09-23

### Fixed
- **Contig Filtering Regex**: Resolved variant calling failure on non-human species by conditioning human-only chromosome filtering (`chr[0-9XYM]+`) to `params.species == 'human'`.
- **Prokaryotic QC Alignment**: Set default `busco_lineage = "bacteria_odb10"` in `conf/prokaryote.config`.
- **JBrowse Subworkflow Coupling**: Decoupled `ANNOTATION_VIS` execution to run strictly when `reference_based` mode with valid `--ref_nuc` is specified.

---

## [0.2.0] - 2026-09-23

### Added
- **nf-core Standardization**:
  - Modularized tool definitions into atomic `modules/nf-core/<tool>/` (`main.nf`, `environment.yml`, `meta.yml`).
  - Structured subworkflows under `subworkflows/local/<subworkflow>/` (`main.nf`, `meta.yml`).
  - Runtime software version collection (`versions.yml`) aggregated into MultiQC.
  - Parameter validation and CLI schema via `nextflow_schema.json`.
  - Centralized tool arguments and publish directives in `conf/modules.config`.
  - Added repository files: `.gitignore`, `.gitattributes`, `.editorconfig`, `.nf-core.yml`, `modules.json`, `tower.yml`, `nf-test.config`, `LICENSE` (MIT), and `CITATIONS.md`.
  - Comprehensive documentation: `docs/usage.md` and `docs/output.md`.
- **Branding**: Rebranded pipeline to **`longvarmeth`** with custom lowercase ASCII art startup banner.

---

## [0.1.2] - 2026-07-13

### Fixed
- **MultiQC File Collision**: Staged NanoPlot, QUAST, and BUSCO output files in sample-specific subdirectories (`nanoplot_${sample_id}`, `quast_${sample_id}`, `busco_${sample_id}`) to prevent Nextflow staging collisions.
- **Sniffles Tabix Indexing**: Standardized bgzip compression output for structural variant VCFs to ensure seamless downstream indexing.
- **HPC Air-Gapped Cache**: Added local data directory detection in SnpEff/VEP download checks to support execution on compute nodes without outbound internet.

---

## [0.1.1] - 2026-07-10

### Added
- **Variant Effect Annotation**: Integrated SnpEff (`snpeff.nf`) and Ensembl VEP (`vep.nf`) downstream of BCFtools small variant calling.
- **Annotation Extraction**: Added `bin/extract_annotations.py` to extract researcher-targeted gene variants (*MEFV*, *NLRP3*, *TNFRSF1A*, *MVK*, *TNFAIP3*).
- **Target Coordinates**: Verified genomic coordinate mapping for *MEFV* familial Mediterranean fever mutations on GRCh38.

---

## [0.1.0] - 2026-07-08

### Added
- **Multi-Species Profiles**: Added predefined profile configurations for human (`conf/human.config`), eukaryotes (`conf/eukaryote.config`), prokaryotes (`conf/prokaryote.config`), and viruses (`conf/virus.config`).
- **Allele-Specific Phasing**: Integrated WhatsHap for variant phasing and read haplotagging (`HP:i:1`, `HP:i:2`).
- **Structural Variants & CNVs**: Integrated Sniffles2 for SV calling and CNVkit for whole-genome copy number alteration profiling.
- **Methylation Profiling**: Integrated Modkit pileup for 5mC/5hmC BedMethyl extraction.
- **Interactive Visualization**: Integrated JBrowse 2 web directory generation and Methylartist locus-specific plotting.
- **Container Version Lock**: Added `conf/version_lock.config` and `versions.yml` to lock exact biocontainer SHA tags.

---

## [0.0.1] - 2026-07-01

### Added
- **Initial Pipeline Architecture**:
  - Nextflow DSL2 core workflow structure.
  - Host GPU-accelerated basecalling and demultiplexing using ONT Dorado.
  - Direct BAM merging and FASTQ extraction with Samtools.
  - Initial read QC with FastQC and Cramino.
  - Initial de novo assembly pipeline with Filtlong, Flye, and Medaka consensus polishing.
  - Reference alignment using Minimap2 and base modification tag transfer via `transfer_mod_tags.py`.
  - Baseline small variant calling with BCFtools mpileup & call.
