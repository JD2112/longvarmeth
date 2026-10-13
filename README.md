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


## Analytical Architecture

### Platform-Specific Workflow Execution Matrix

| Stage | Oxford Nanopore <br> (`--platform ont`) | PacBio HiFi <br> (`--platform pacbio`) | Inputs $\rightarrow$ Outputs |
| :--- | :--- | :--- | :--- |
| **Input Checksum** | Ingestion SHA-256 verification | Ingestion SHA-256 verification | Raw files $\rightarrow$<br>Audit manifest |
| **Read QC** | `FastQC` + `NanoPlot`<br>`Cramino` + `samtools` | `Cramino` + `FastQC`<br>`NanoPlot` | Reads / BAM $\rightarrow$<br>Yield & Q-scores |
| **Basecalling**<br>*(ONT raw)* | Host GPU `Dorado`<br>+ `pod5 convert` demux | *On-instrument demux*<br>(Revio / Sequel IIe) | `pod5` / `fast5` $\rightarrow$<br>Basecalled BAM |
| **Assembly**<br>*(de novo)* | `Flye --nano-hq`<br>+ `Medaka` polishing | `Hifiasm` (diploid)<br>+ `Racon` polishing | Long reads $\rightarrow$<br>Assembly FASTA + BUSCO |
| **Alignment** | `minimap2 -ax map-ont`<br>+ `MM`/`ML` tag restoration | `pbmm2 align --preset HIFI`<br>*(native kinetic preservation)* | Reads + Ref $\rightarrow$<br>Sorted BAM (`MM`/`ML`) |
| **Small Variants** | `Clair3` (ONT model)<br>`DeepVariant` / `BCFtools` | `DeepVariant` (HiFi model)<br>`Clair3` / `BCFtools` | Sorted BAM $\rightarrow$<br>SNV / InDel VCF |
| **SVs, CNVs & STRs** | `Sniffles2` (SVs)<br>`CNVkit` (CNVs) | `Sniffles2` (SVs) + `HiFiCNV`<br>`TRGT` (tandem repeats) | Sorted BAM $\rightarrow$<br>SV, CNV & STR calls |
| **Epigenomics** | `modkit pileup`<br>(5mC / 5hmC / 6mA motifs) | `pbcpgtools` + `methbat`<br>+ `modkit pileup` | Sorted BAM $\rightarrow$<br>BedMethyl & CpG scores |
| **Phasing** | `WhatsHap phase`<br>+ `haplotag` (`HP:i:1/2`) | `HiPhase`<br>*(co-phases SNV + SV + STR)* | VCFs + BAM $\rightarrow$<br>Phased VCF & Haplotagged BAM |
| **Annotation** | `Ensembl VEP`<br>`SnpEff` + ClinVar | `Ensembl VEP`<br>`SnpEff` + ClinVar | Phased VCF $\rightarrow$<br>Annotated clinical VCF |
| **Governance** | `SANITIZE_METADATA` (PHI scrub)<br>+ cryptographic SHA-256 manifest | `SANITIZE_METADATA` (PHI scrub)<br>+ cryptographic SHA-256 manifest | Artifacts $\rightarrow$<br>Sanitized BAM & manifest |
| **Visualization** | `JBrowse 2` web workspace<br>`Methylartist` locus plots<br>Interactive Methylation Dashboard | `JBrowse 2` web workspace<br>`Methylartist` locus plots<br>Interactive Methylation Dashboard | BAMs, VCFs $\rightarrow$<br>Web workspace & HTML reports |
| **Reporting** | `MultiQC` consolidated dashboard<br>+ Executive clinical summary | `MultiQC` consolidated dashboard<br>+ Executive clinical summary | Metrics $\rightarrow$<br>MultiQC HTML report |

### Multi-Organism Tailoring Matrix

| Organism Profile | Target Scope | Ploidy | Specialized Calling & Assembly Behavior |
| :--- | :--- | :---: | :--- |
| **Human**<br>(`-profile human`) | Clinical diagnostics, rare disease,<br>germline/somatic WGS | Diploid<br>(`2`) | • Full variant matrix (`DeepVariant` + `Sniffles2` + `TRGT` + `HiPhase`)<br>• Whole-genome 5mCG methylation & imprinting DMRs<br>• Clinical VEP/ClinVar consequence annotation & PHI scrubbing |
| **Eukaryote**<br>(`-profile eukaryote`) | Non-model plants, animals,<br>fungi, invertebrates | Diploid<br>(`2`) | • Reference-guided or *de novo* assembly with `Hifiasm` / `Flye`<br>• Lineage-specific gene completeness benchmark via `BUSCO`<br>• Reference-free phased methylomics |
| **Prokaryote**<br>(`-profile prokaryote`) | Bacteria, archaea, plasmids,<br>microbial isolates | Haploid<br>(`1`) | • Single-contig circularized assembly via `Flye` / `Hifiasm`<br>• Circularity verification & replication origin analysis<br>• Restriction-modification (R-M) motif detection (4mC, 6mA, 5mC) |
| **Virus**<br>(`-profile virus`) | Viral genomics, amplicons,<br>quasispecies, phages | Haploid<br>(`1`) | • Amplicon-aware quality & length filtering<br>• Low-frequency intra-host diversity variant detection (relaxed 5% threshold)<br>• Metagenomic assembly mode (`--meta`) |

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
| `--input` | string | `samplesheet.csv` | Path to CSV samplesheet (`run_id,input_path,input_type,...`). Supports `pod5`, `fast5`, `ubam`, `pacbio_bam`, or pre-aligned `aligned_bam`. |
| `--outdir` | string | `results` | Output directory destination |
| `--platform` | string | `ont` | Long-read sequencing technology (`ont`, `pacbio`) |
| `--species` | string | `eukaryote` | Target biology profile (`human`, `eukaryote`, `prokaryote`, `virus`) |
| `--run_type` | string | `reference_based` | Analysis workflow mode (`reference_based`, `de_novo_assembly`) |
| `--ref_nuc` | string | `""` | Path to nuclear reference genome FASTA |
| `--snv_caller` | string | `bcftools` | Small variant caller (`bcftools`, `deepvariant`, `clair3`) |
| `--sv_caller` | string | `sniffles` | Structural variant caller (`sniffles`, `severus`, `pbsv`) |
| `--phasing_tool` | string | `whatshap` | Haplotype phasing engine (`whatshap`, `hiphase`) |
| `--cnv_caller` | string | `cnvkit` | Copy number variant caller (`cnvkit`, `hificnv`) |
| `--skip_variants` | boolean | `false` | Skip all variant calling (SNVs, SVs, CNVs, STRs) to run a pure methylation pipeline |
| `--skip_methylation` | boolean | `false` | Skip all epigenomics (Modkit / pb-cpg-tools) to run a pure variant pipeline |
| `--skip_alignment` | boolean | `false` (auto for `aligned_bam`) | Bypass Minimap2/pbmm2 alignment when providing pre-aligned coordinate-sorted BAMs |
| `--enable_trgt` | boolean | `false` | Enable tandem repeat genotyping via TRGT |
| `--repeat_bed` | string | `""` | Path to TRGT repeat definition BED catalog |
| `--enable_phi_sanitization`| boolean | `true` | Privacy-by-default: strip personal health information (PHI) & hospital run barcodes from BAM/VCF headers |
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
