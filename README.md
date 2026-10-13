# longvarmeth

[![quindecagon](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/quindecagon.json)](https://github.com/JD2112/quindecagon)
[![nf-core lint](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/nfcore_lint.json)](https://nf-co.re)
[![SBOM](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/syft.json)](https://github.com/anchore/syft)
[![secrets](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/secrets.json)](https://github.com/gitleaks/gitleaks)
[![code style](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/black.json)](https://github.com/psf/black)
[![flake8](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/flake8.json)](https://flake8.pycqa.org)
[![R audit](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JD2112/longvarmeth/badges/r_audit.json)](https://github.com/r-lib/lintr)

A modular, clinical-grade Nextflow DSL2 platform for single-molecule long-read genomics and epigenomics (**Oxford Nanopore & PacBio HiFi**). It features an explicit 12-stage analytical architecture with host GPU-accelerated basecalling/demux, pre- and post-alignment QC, *de novo* assembly & polishing, kinetic tag preservation (`MM`/`ML`), small and structural variant calling (SNV/InDel/SV/CNV/STR), haplotype co-phasing, strict clinical governance & PHI sanitization (inspired by `milou`), and interactive reporting.

## Metro Architecture

![longvarmeth Metro Diagram](./docs/pipeline.png)


## 12-Stage Analytical Architecture

### Platform-Specific Workflow Execution Matrix

| Stage | Oxford Nanopore (`--platform ont`) | PacBio HiFi (`--platform pacbio`) | Inputs $\rightarrow$ Outputs |
| :--- | :--- | :--- | :--- |
| **1. Input Checksum** | Ingestion SHA-256 verification | Ingestion SHA-256 verification | Raw files $\rightarrow$ Audit manifest |
| **2. Read QC** | `FastQC` + `Cramino` + `NanoPlot` | `Cramino` + `FastQC` + `NanoPlot` | BAMs / FASTQ $\rightarrow$ Yield & Q-scores |
| **3. Basecalling** *(ONT raw)* | Host GPU `Dorado` basecalling & demux | *On-instrument demux (Revio / Sequel)* | `pod5`/`fast5` $\rightarrow$ Basecalled BAM |
| **3. Assembly** *(de novo)* | `Flye --nano-hq` + `Medaka` polishing | `Hifiasm` diploid + `Racon` polishing | Reads $\rightarrow$ Polished FASTA + BUSCO |
| **4. Alignment** | `minimap2 -ax map-ont` + tag restoration | `pbmm2 align --preset HIFI` (native) | Reads + Ref $\rightarrow$ Sorted BAM (`MM`/`ML`) |
| **5. Small Variants** | `Clair3` (ONT model) / `DeepVariant` / `BCFtools` | `DeepVariant` (HiFi model) / `Clair3` | Sorted BAM $\rightarrow$ SNV/InDel VCF |
| **6. SVs, CNVs & STRs** | `Sniffles2` (SVs) + `CNVkit` (CNVs) | `Sniffles2` + `HiFiCNV` (CNVs) + `TRGT` (STRs)| Sorted BAM $\rightarrow$ SV, CNV, and STR calls |
| **7. Epigenomics** | `modkit pileup` (5mC/5hmC + 6mA motifs) | `pbcpgtools` + `methbat` + `modkit pileup` | Sorted BAM $\rightarrow$ BedMethyl & CpG scores |
| **8. Phasing** | `WhatsHap phase` & `haplotag` (`HP:i:1/2`) | `HiPhase` (co-phases SNV + SV + STR) | VCFs + BAM $\rightarrow$ Phased VCF & BAM |
| **9. Annotation** | `Ensembl VEP` + `SnpEff` + ClinVar | `Ensembl VEP` + `SnpEff` + ClinVar | Phased VCF $\rightarrow$ Annotated VCF |
| **10. Governance** | PHI / hospital barcode sanitization & SHA-256 | PHI / hospital barcode sanitization & SHA-256 | Artifacts $\rightarrow$ Sanitized BAM + manifest |
| **11. Visualization**| `JBrowse 2` tracks, `Methylartist` & **Interactive Methylation Dashboard** | `JBrowse 2` tracks, `Methylartist` & **Interactive Methylation Dashboard** | BAMs, VCFs $\rightarrow$ Web workspace & HTML reports |
| **12. Reporting** | `MultiQC` + Executive Clinical Summary | `MultiQC` + Executive Clinical Summary | Metrics $\rightarrow$ MultiQC HTML report |

### Multi-Organism Tailoring Matrix

| Organism Profile | Target Scope | Ploidy | Specialized Calling & Assembly Behavior |
| :--- | :--- | :---: | :--- |
| **Human** (`-profile human`) | Clinical diagnostics, rare disease, germline/somatic WGS | Diploid (`2`) | Full variant matrix (`DeepVariant` + `Sniffles2` + `TRGT` + `HiPhase`), whole-genome 5mCG methylation & imprinting DMRs, clinical VEP/ClinVar annotation. |
| **Eukaryote** (`-profile eukaryote`) | Non-model plants, animals, fungi, invertebrates | Diploid (`2`) | Reference-guided or *de novo* assembly with `Hifiasm`/`Flye`, lineage-specific gene completeness benchmark via `BUSCO`, and reference-free phased methylomics. |
| **Prokaryote** (`-profile prokaryote`) | Bacteria, archaea, plasmids, microbial isolates | Haploid (`1`) | Single-contig circularized genome assembly via `Flye`/`Hifiasm`, circularity verification, and bacterial restriction-modification (R-M) motif detection (4mC, 6mA, 5mC). |
| **Virus** (`-profile virus`) | Viral genomics, amplicons, quasispecies, phages | Haploid (`1`) | Amplicon-aware quality filtering, low-frequency intra-host diversity variant detection (relaxed 5% threshold), and metagenomic assembly (`--meta`). |

## Quickstart

### 1. Oxford Nanopore Reference-Based Run
```bash
nextflow run main.nf \
  -profile human,singularity \
  --platform ont \
  --input samplesheet.csv \
  --ref_nuc /path/to/GRCh38.fa \
  --outdir results_ont
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
  --outdir results_pacbio
```

### 3. Bacterial Assembly & Methylation
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



## Key Parameters

You can configure execution parameters from the command line or via custom config files:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input` | string | `samplesheet.csv` | Path to CSV samplesheet (`run_id,input_path,input_type,kit,model,modified_bases_model`) |
| `--outdir` | string | `results` | Output directory destination |
| `--platform` | string | `ont` | Long-read sequencing technology (`ont`, `pacbio`) |
| `--species` | string | `eukaryote` | Target biology profile (`human`, `eukaryote`, `prokaryote`, `virus`) |
| `--run_type` | string | `reference_based` | Analysis workflow mode (`reference_based`, `de_novo_assembly`) |
| `--ref_nuc` | string | `""` | Path to nuclear reference genome FASTA |
| `--snv_caller` | string | `bcftools` | Small variant caller (`bcftools`, `deepvariant`, `clair3`) |
| `--sv_caller` | string | `sniffles` | Structural variant caller (`sniffles`, `severus`, `pbsv`) |
| `--phasing_tool` | string | `whatshap` | Haplotype phasing engine (`whatshap`, `hiphase`) |
| `--cnv_caller` | string | `cnvkit` | Copy number variant caller (`cnvkit`, `hificnv`) |
| `--enable_trgt` | boolean | `false` | Enable tandem repeat genotyping via TRGT |
| `--repeat_bed` | string | `""` | Path to TRGT repeat definition BED catalog |
| `--enable_phi_sanitization`| boolean | `false` | Strip personal health information (PHI) & hospital run IDs from BAM/VCF headers |
| `--enable_audit_trail` | boolean | `true` | Generate cryptographic SHA-256 manifests and forensic provenance records |
| `--enable_vep` | boolean | `false` | Enable Ensembl VEP functional variant annotation |
| `--enable_snpeff` | boolean | `false` | Enable SnpEff variant annotation |
| `--enable_methylation_dashboard` | boolean | `false` (auto in `human`) | Generate interactive AG-Grid & IGV-Web clinical methylation report (adapted from `gms-nallo`) |
| `--dorado_bin` | string | `/path/to/dorado`| Host path to Dorado binary for host GPU-accelerated basecalling |



## Species Profiles

Preconfigured profiles load ploidy, reference genomes, and calling parameters:

| Profile | Target Organism Type | Key Settings |
| :--- | :--- | :--- |
| `-profile human` | *Homo sapiens* | Diploid (`--ploidy 2`), hg38 defaults, WhatsHap/HiPhase phasing, VEP/SnpEff enabled |
| `-profile eukaryote` | Non-Model Eukaryotes | Diploid ploidy, reference-guided or Hifiasm/Flye assembly with BUSCO lineage QC |
| `-profile prokaryote`| Bacteria & Archaea | Haploid (`--ploidy 1`), circularized assembly, bacterial 4mC/6mA/5mC motif detection |
| `-profile virus` | Viruses (Amplicons / Quasispecies) | Haploid (`--ploidy 1`), relaxed 5% frequency cutoffs for intra-host diversity |



## Testing

Run test suites using official, lightweight public `nf-core/test-datasets`:

```bash
# Nanopore test run:
nextflow run main.nf -profile test_ont,singularity

# PacBio HiFi test run:
nextflow run main.nf -profile test_pacbio,singularity

# Dry-run syntax & topology check:
nextflow run main.nf -profile test_pacbio -stub
```



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

## Reproducibility & Version Lock

To ensure exact forensic reproducibility across clinical production runs, the pipeline includes a **`version_lock`** profile. When activated, it forces Nextflow to use exact container SHA digests listed in `conf/version_lock.config` and documented in `versions.yml`.

```bash
nextflow run main.nf -profile human,singularity,version_lock [parameters]
```

## Citations

If you use `longvarmeth` in your research, please cite the pipeline:

> Das, J. (2026). *longvarmeth: A Unified Nextflow DSL2 Platform for Clinical Long-Read Genomics and Epigenomics (Nanopore & PacBio HiFi)*. GitHub: [https://github.com/JD2112/longvarmeth](https://github.com/JD2112/longvarmeth).

Please also cite the underlying bioinformatics tools listed in [CITATIONS.md](CITATIONS.md).


## License & Acknowledgements

* **License**: This project is licensed under the [MIT License](LICENSE).
* **Acknowledgements**:
  * Built using community best practices from [nf-core](https://nf-co.re).
  * Clinical governance and audit trail architecture adapted from [`milou`](https://github.com/JD2112/milou).
  * Interactive clinical methylation dashboard adapted from [`gms-nallo-methylation-dashboard`](https://github.com/JD2112/gms-nallo-methylation-dashboard).
  * Metro-map architecture rendered with [nf-metro](https://github.com/seqeralabs/nf-metro).
