# FastqToGeneCounts

[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

A Snakemake workflow for processing bulk RNA-seq data from NCBI GEO accession numbers or FASTQ files.
This pipeline performs read quality control, alignment to a reference genome, and generates gene count matrices suitable
for downstream analysis such as differential expression or metabolic modeling.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Conda Environment Setup](#conda-environment-setup)
- [Configuration](#configuration)
    - [Sample File (MASTER_CONTROL)](#sample-file-master_control)
    - [Output Settings](#output-settings)
    - [Processing Options](#processing-options)
    - [Genome Settings](#genome-settings)
- [Profile Configuration](#profile-configuration)
    - [SLURM Cluster](#slurm-cluster)
    - [Local Execution](#local-execution)
- [Running the Pipeline](#running-the-pipeline)
    - [Dry Run](#dry-run)
    - [Full Execution](#full-execution)
- [Output Description](#output-description)
- [Benchmarking](#benchmarking)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

## Overview

This workflow:

1. Downloads SRA files from NCBI/GEO (if using prefetch) or reads local FASTQ files
2. Converts SRA files to FASTQ format using parallel-fastq-dump
3. Performs quality control with FastQC before and after trimming
4. (Optionally) Trims adapters and low-quality bases with Trim Galore
5. Aligns reads to a chosen reference genome using STAR
6. Quantifies gene counts using Salmon
7. Collects RNA-seq metrics with Picard
8. Calculates insert sizes and fragment sizes
9. Screens for common contaminants
10. Generates a MultiQC report combining all quality metrics

## Installation

Clone the repository using Git:

```bash
git clone https://git.unl.edu/HelikarLab/FastqToGeneCounts.git
cd FastqToGeneCounts
```

## Conda Environment Setup

Create and activate a new Conda environment with the required dependencies:

```bash
conda create -n fastq2genecounts python=3.10 snakemake
conda activate fastq2genecounts

# If you are running this pipeline and want to use SLURM, install the executor plugin:
pip install snakemake-executor-plugin-slurm
```

Verify the installation:

```bash
snakemake --version
```

## Setup

### Configuration File

The `config.yaml` file contains all customization options available to the pipeline. The options in this file are
explained below.

### Sample File (`MASTER_CONTROL`)

The `MASTER_CONTROL` setting specifies the CSV file containing sample information. This file must have the
following four columns as a header row:

| Column        | Description         | Valid Values                                                                                                                                |
|---------------|---------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `srr`         | SRA accession code  | SRR###### (leave this empty if you are using local files)                                                                                   |
| `sample`      | Sample identifier   | Prefix with any alphanumeric string (e.g., tissue name), plus the "Study", "Run", and "Replicate" for that sample, separated by underscores |
| `endtype`     | Sequencing type     | `PE` (paired-end) or `SE` (single-end)                                                                                                      |
| `prep_method` | Library preparation | `total` (total RNA) or `mrna` (PolyA/mRNA)                                                                                                  |

Example CSV file:

```csv
srr,sample,endtype,prep_method
SRR101,tissueA_S1R1,PE,mrna
SRR102,tissueA_S1R2,SE,total
SRR103,tissueB_S1R1r1,PE,total
SRR104,tissueB_S1R1r2,PE,total
```

For local FASTQ files (see [Processing Options](#processing-options) below), leave the srr column empty. If the
sampe uses paired-end sequencing, the pipeline will automatically find the forward and reverse reads based on the
sample name provded in the file (e.g., `tissueA_S1R1_1.fastq.gz` and `tissueA_S1R1_2.fastq.gz`)

```csv
srr,sample,endtype,prep_method
,tissueA_S1R1,PE,total
,tissueA_S1R2,SE,total
```

### Output Settings

| Setting           | Description                                                                                    |
|-------------------|------------------------------------------------------------------------------------------------|
| `ROOTDIR`         | Directory where results will be saved (default: `results`)                                     |
| `EXPERIMENT_NAME` | Optional subdirectory name for organizing experiment outputs                                   |
| `BENCHMARK_DIR`   | Directory for benchmark output (default: `benchmarks`)                                         |
| `BENCHMARK_TIMES` | Number of times to run each rule for benchmarking (default: 1, minimal benchmarking performed) |
| `LOG_DIR`         | Directory for log files (default: `logs`)                                                      |

### Processing Options

| Setting                       | Description                                  | Default |
|-------------------------------|----------------------------------------------|---------|
| `PERFORM_DUMP_FASTQ`          | Download + convert SRA files to FASTQ format | `True`  |
| `PERFORM_TRIM`                | Trim adapters and low-quality bases          | `True`  |
| `PERFORM_SCREEN`              | Screen for contamination                     | `True`  |
| `PERFORM_GET_RNASEQ_METRICS`  | Collect RNA-seq metrics with Picard          | `True`  |
| `PERFORM_GET_INSERT_SIZE`     | Calculate insert sizes                       | `True`  |
| `PERFORM_GET_FRAGMENT_SIZE`   | Calculate fragment sizes with RSeQC          | `True`  |
| `BYPASS_REPLICATE_VALIDATION` | Skip replicate count validation (for COMO)   | `False` |
| `BYPASS_GENOME_VALIDATION`    | Skip genome file validation before running   | `False` |

#### Local FASTQ Mode

To process local FASTQ files instead of downloading from SRA:

1. Set `PERFORM_DUMP_FASTQ: False`
2. Set `LOCAL_FASTQ_FILES: "<dirname>"` (change `<dirname>` to the directory path)
3. Place your FASTQ files in the specified directory
4. Leave the `srr` column empty in your MASTER_CONTROL CSV

### Genome Settings

The pipeline automatically downloads and prepares reference genome files from Ensembl. Configure these settings under
the `GENOME` section:

| Setting         | Description                | Example                          |
|-----------------|----------------------------|----------------------------------|
| `SAVE_DIR`      | Directory for genome files | `genome`                         |
| `TAXONOMY_ID`   | NCBI taxonomy ID           | `9606` (human), `10090` (mouse)  |
| `VERSION`       | Ensembl release version    | `115`, `latest`, `release-112`   |
| `SHOW_PROGRESS` | Show download progress     | `True`                           |
| `FASTA_VERSION` | FASTA file type            | `primary_assembly` or `toplevel` |

Find your taxonomy ID at https://www.ncbi.nlm.nih.gov/Taxonomy

Available Ensembl releases are listed at https://www.ftp.ensembl.org/pub/

## Profile Configuration

Snakemake profiles configure how jobs are submitted to your compute environment. Two profiles are included with this
pipeline.

### SLURM Cluster

Use the SLURM profile for high-performance computing clusters:

```bash
snakemake --profile profiles/cluster
```

The SLURM profile configuration (`profiles/cluster/config.v8+.yaml`) includes:

| Setting         | Description                                | Default |
|-----------------|--------------------------------------------|---------|
| `cores`         | Maximum cores used in workflow             | 10      |
| `jobs`          | Maximum concurrent jobs to submit to SLURM | 250     |
| `restart-times` | Retry attempts on failure                  | 0       |
| `use-conda`     | Enable Conda environments                  | `True`  |

> [!WARNING] Conda
> The `use-conda` setting should always be kept to **True**, otherwise Snakemake assumes all required dependencies
> are installed in the current environment. By keeping this value as True, Snakemake will create the required conda
> environments and activate them for each required rule.

#### Customizing SLURM Account

Edit `profiles/cluster/config.v8+.yaml` to change the SLURM account:

```yaml
sbatch
  --job-name=smk-{rule}-{wildcards}
  --account=YOUR_ACCOUNT_NAME  # Change to your account name: `sacctmgr show user $USER accounts`
  --cpus-per-task={threads}
...
```

### Local Execution

For execution on a single machine without a cluster scheduler:

```bash
snakemake --profile profiles/local
```

> [!warning]
> Local execution may be slow for large datasets and, unless you are running this pipeline on a workstation, is not 
> recommended for genome alignment due to memory constraints.

### Additional Parameters
If you would like to change any of the default parameters in the profile, you can either:
1. Modify the profile to your liking; this will become the new default for all future pipeline workflows, or
2. Override the parameters on a per-run basis by including them on the command line. Execute `snakemake --help` to 
   view all options

## Running the Pipeline

### Dry Run

A dry run can be performed to verify configuration:

```bash
snakemake --profile profiles/cluster --dry-run # or `--profile profiles/local` if using a local workstation
```

The dry run will:

- Validate the Snakefile syntax
- Print the defined output directories
- Display all rules that will be executed
- Check that configuration is properly set up
- Verify input files are accessible

### Full Execution

After confirming the dry run succeeds, run the full pipeline:

```bash
snakemake --profile profiles/cluster  # or `--profile profiles/local` if using a local workstation
```

#### Keeping the Session Alive

Since Snakemake does not daemonize by default, use `screen` or `tmux` before running Snakemake to keep the 
main process running:

- [Screen](https://www.gnu.org/software/screen/manual/screen.html)
- [Tmux](https://man7.org/linux/man-pages/man1/tmux.1.html)

> [!NOTE] Screen vs Tmux
> Screen is easer to set up and use, but Tmux offers more powerful features

```bash
# Start a screen session
screen -S snakemake

# Detach from session: Ctrl+a, d

# Reattach to session
screen -r snakemake
```

#### Log Files

Log files are stored in the `logs` directory organized by tissue and rule:
```
logs/  # or the defined `LOG_DIR` configuration option
  {tissue}/
    {rule}/
      {tissue}_{rule}.log
```

## Benchmarking

The pipeline supports benchmarking to measure execution time for each rule. To enable benchmarking:

1. Set `BENCHMARK_TIMES` to your desired number of repetitions (e.g., `3`)
2. Run the pipeline normally

Generate benchmark plots using the provided script:

```bash
python utils/benchmark.py
```

This script reads benchmark data from the `benchmarks/` directory and produces visualization plots.

> [!WARNING] Bencharmking
> The pipeline will run `BENCHMARK_TIMES` *for each input file*. If set to 2, the entire pipeline will run twice to
> collect benchmark statistics.

## Getting Help

Please [open an issue](https://git.unl.edu/helikarlab/FastqToGeneCounts/issues) with:

- The exact error message
- Your configuration file
- The rule and sample that failed
- Relevant log files from `logs/`

## Citation

If you use this pipeline in your research, please cite it as:

```bibtex
@software{FastqToGeneCounts,
  author = {Josh Loecker, Brandt Bessell, Bhanwar lal Puniya, Tomáš Helikar},
  title = {FastqToGeneCounts: Automated Bulk RNA-seq Alignment},
  url = {https://git.unl.edu/helikarlab/FastqToGeneCounts}
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a merge request for:

- Bug reports
- Feature requests
- Documentation improvements
- Code enhancements
