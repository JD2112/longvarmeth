# longvarmeth

A modular, robust Nextflow DSL2 pipeline for long-read (ONT) sequencing data. It supports GPU-accelerated basecalling/demux on the host, QC, *de novo* assembly, polishing, reference alignment, methylation tag transfer, variant calling (SNV/InDel/SV/CNV), allele-specific phasing, and interactive visualization (JBrowse 2 / Methylartist / MultiQC).

## Key Features
* **Basecalling & Demultiplexing**: Host GPU-accelerated basecalling and barcode demultiplexing using ONT **Dorado**.
* **Quality Control & Polishing**: Multi-tool read quality assessment (**FastQC**, **Cramino**, **NanoPlot**), *de novo* assembly (**Flye**), assembly evaluation (**QUAST**, **BUSCO**), and consensus polishing (**Medaka**).
* **Alignment & Modification Transfer**: Reference-based alignment with **Minimap2** and base modification tag preservation (`MM`/`ML` tags).
* **Variant Calling & Annotation**:
  * **SNV / InDel**: Small variant discovery via **BCFtools**.
  * **Structural Variants**: High-sensitivity SV calling with **Sniffles2**.
  * **Copy Number Variants**: Whole-genome copy-number profiling via **CNVkit**.
  * **Functional Annotation**: Variant effect prediction with **SnpEff** and Ensembl **VEP**.
* **Phasing & Methylation Profiling**:
  * Read haplotagging and variant phasing via **WhatsHap**.
  * 5mC / 5hmC modification pileup extraction using ONT **Modkit**.
* **Interactive Visualization**: Automatically generates **JBrowse 2** interactive genome browser tracks and locus-specific modification plots with **Methylartist**.
* **Multi-Species Support**: Pre-configured species profile flags (`-profile human`, `-profile eukaryote`, `-profile prokaryote`, `-profile virus`).

## Quickstart

Run the pipeline using a specific species and container engine profile:

* **Human Genome WGS/WES (reproducible version locked)**:
```bash
nextflow run main.nf \
  -profile human,singularity,version_lock \
  --enable_snpeff true \
  --snpeff_cache /path/to/snpEff
```

* **Bacterial / Prokaryotic De Novo Assembly (Flye + Medaka + QUAST + BUSCO)**:
```bash
nextflow run main.nf \
  -profile prokaryote,singularity \
  --input samplesheet.csv \
  --run_type de_novo_assembly
```

* **Bacterial Reference-Based Alignment & Variant / Methylation Calling**:
```bash
nextflow run main.nf \
  -profile prokaryote,singularity \
  --input samplesheet.csv \
  --run_type reference_based \
  --ref_nuc /path/to/bacteria.fasta
```

* **Other Eukaryotes (e.g., Fungi, Plants, Invertebrates)**:
```bash
nextflow run main.nf \
  -profile eukaryote,singularity \
  --input samplesheet.csv \
  --ref_nuc /path/to/genome.fasta \
  --busco_lineage fungi_odb10
```

## Species Profiles

Species profiles load predefined reference genomes, ploidy, and filtering criteria. Select one species profile for your run:

| Profile | Target Organism Type | Key Settings |
| :--- | :--- | :--- |
| `-profile eukaryote` | Eukaryotes (Artemia sp. defaults) | Diploid ploidy, default Artemia nuclear & mtDNA references |
| `-profile human` | Homo sapiens | hg38 reference (`/data/references/hg38.fa`), diploid ploidy, WhatsHap phasing enabled |
| `-profile prokaryote`| Bacteria & Archaea | Haploid ploidy, single circular reference support |
| `-profile virus` | Viruses | Haploid ploidy, small circular/linear genome support |

## Engine Profiles

Configure how pipeline processes run their official biocontainers:

* **Singularity**: `-profile singularity` (Recommended for HPC, includes `-B /data:/data` bind mount option)
* **Apptainer**: `-profile apptainer` (Modern HPC systems)
* **Docker**: `-profile docker` (Local developer workstations)
* **Conda**: `-profile conda` (Fallback local environment activation)

## Reproducibility & Version Lock

To ensure exact process reproducibility across runs, the pipeline includes a **`version_lock`** profile. When activated, it forces Nextflow to use exact container hashes/tags listed in `conf/version_lock.config` and documented in `versions.yml`.

To run with the version lock active:
```bash
nextflow run main.nf -profile human,singularity,version_lock [parameters]
```

## Key Parameters

You can override any parameter defined in `conf/params.config` directly from the command line:

### Input / Output Options
* `--input <path>`: Path to input samplesheet CSV (default: `samplesheet.csv`).
* `--outdir <path>`: Directory to write output files (default: `results`).

### Run Modes
* `--run_type <mode>`:
  * `reference_based` (default): Align reads directly to references and call variants/methylation.
  * `de_novo_assembly`: Assemble reads with *Flye*, polish with *Medaka*, and run *QUAST*/*BUSCO* assembly QC.

### Feature Toggles
* `--enable_phasing <true|false>`: Enable whatsHap read phasing and haplotagging (default: `false` but defaults to `true` for `human`).
* `--methylartist_locus <coordinates>`: Genomic coordinates to plot locus-specific methylation with *Methylartist* (e.g. `--methylartist_locus "NC_001807:1-1000"`).

### Dorado GPU Basecalling Options
* `--dorado_bin <path>`: Path to local dorado binary on the host.
* `--gpu_devices <device>`: Dorado CUDA device select flag (e.g. `--gpu_devices "cuda:all"`).

## Verification / Dry Run

To check nextflow syntax and dry-run/stub process steps without loading large files:

```bash
nextflow run main.nf -c tests/test.config -stub
```

To view JBrowse, the directory must be served through a local or remote web server. 

Here are the two easiest ways to view it:

### Option A: If you downloaded the folder to your local computer
```bash
#1. Open your terminal, and navigate (`cd`) into the JBrowse sample directory:
cd /path/to/results/jbrowse/sample_bc05
#2. Start a simple Python web server:
python3 -m http.server 8000
python3 -c "import http.server, mimetypes; mimetypes.init(); mimetypes.encodings_map.pop('.gz', None); http.server.test(HandlerClass=http.server.SimpleHTTPRequestHandler, port=8000)"
#3. Open your browser and navigate to: **`http://localhost:8000`**
```

### Option B: If the files are on your remote server (via SSH)
```bash
#1. On the remote server, navigate to the directory and start the server:

cd results/jbrowse/sample_bc05
python3 -m http.server 8000

#2. On your local machine, open a new terminal and set up an SSH port-forwarding tunnel:
ssh -L 8000:localhost:8000 username@your-server.domain.com

#3. Open your browser and navigate to: **`http://localhost:8000`**
```

### Option A: If `bgzip` and `tabix` are installed on your server path
```bash
# Navigate to the folder containing the sniffles SV vcf file
cd results/variants/svs/

# Decompress standard gzip file
gunzip sample_bc05_svs.vcf.gz

# Re-compress using bgzip
bgzip sample_bc05_svs.vcf

# Index the new bgzipped file (produces .tbi)
tabix -p vcf sample_bc05_svs.vcf.gz
```

### Option B: If you do not have `bgzip` installed on your host system
You can use `apptainer`/`singularity` (which is already configured on your system) to run them from a container:

```bash
cd results/variants/svs/

# Decompress standard gzip
gunzip sample_bc05_svs.vcf.gz

# Run bgzip & tabix via Singularity
singularity exec docker://quay.io/biocontainers/tabix:1.11--h5e330a3_1 bgzip sample_bc05_svs.vcf
singularity exec docker://quay.io/biocontainers/tabix:1.11--h5e330a3_1 tabix -p vcf sample_bc05_svs.vcf.gz
```

---

## Documentation

* [Usage Documentation](docs/usage.md): Detailed parameter descriptions and usage examples.
* [Output Documentation](docs/output.md): Directory layout and descriptions of output files.
* [Citations](CITATIONS.md): Literature citations for pipeline modules and underlying tools.
* [License](LICENSE): MIT License.

