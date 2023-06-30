---
title: Setting up FastqToGeneCounts
sidebar: sidebar
permalink: fastq_setup.html
summary: This is an overview of how to set up the pipeline
last_updated: Sept 22, 2022
---

## Overview
Following the overview is a step-by-step guide on how to set up the pipeline. A brief description is as follows:

1. Create a conda environment containing:
    1. Mamba
    2. Snakemake ([version 7.19.1](https://snakemake.readthedocs.io/en/stable/project_info/history.html#id2))
2. Create a CookieCutter template/profile based on the [Slurm profile](https://github.com/Snakemake-Profiles/slurm)
3. Modification of configuration variables

## Conda Environment
It should be noted that some of these steps can take a considerable length of time, especially depending on your hardware.

Choose a name for your conda environment. This tutorial will be using `snakemake` for its environment name




### 1. Create the Environment
```bash
conda create -n snakemake
```

### 2. Activate the Environment
```bash
conda activate snakemake
```

### 3. Install Mamba
{% include note.html content="Mamba is [recommended by Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as it is much faster than Conda. Additionally, Conda can have issues installing Snakemake dependencies." %}
```bash
conda install -n snakemake -c conda-forge mamba
```

|    Component     |                          Description                          |
|:----------------:|:-------------------------------------------------------------:|
|     `conda`      |                   The command to run Conda                    |
|    `install`     |               The command to install a package                |
|  `-n snakemake`  | The name of the conda environment to install into (snakemake) |
| `-c conda-forge` |           The channel to install from (Conda Forge)           |
|     `mamba`      |                    The package to install                     |


### 4. Install Snakemake
After installing Mamba, we can install Snakemake. This will take a while, as Snakemake has a lot of dependencies. We have specified version `6.15.5`, as version 7 includes breaking changes not yet resolved.
```bash
mamba install -n snakemake -c conda-forge -c bioconda snakemake=7.19.1 pydantic=1.10.3
```

|          Component           |                          Description                          |
|:----------------------------:|:-------------------------------------------------------------:|
|           `mamba`            |                   The command to run Mamba                    |
|          `install`           |               The command to install a package                |
|        `-n snakemake`        | The name of the conda environment to install into (snakemake) |
| `-c conda-forge -c bioconda` |    The channels to install from (Conda Forge and Bioconda)    |
|      `snakemake=6.15.5`      |                    The package to install                     |
|       `pydantic=10.3`        |             Install pydantic for dataclass usage              |


### 5. Test the Snakemake installation
You should receive no errors at this point, and the output should follow the following example.<br>
If you encounter errors, please [open an issue](https://github.com/HelikarLab/FastqToGeneCounts/issues)
```bash
snakemake --version
# Returns `7.25.0`
```

### Profile Benefits
Snakemake Profiles provide multiple benefits:
1. No need to write Slurm scripts
2. Jobs are automatically submitted to the scheduler, and killed if they exceed our resources
3. Snakemake will automatically retry jobs that fail
4. Snakemake will schedule jobs independently, so that they can run in parallel

Perhaps the most important part of Profiles is Point 4. Instead of requesting, for example, 40 cores, 50GB of RAM, and multiple hours for the entire Snakemake workflow (as this is the maximum resources we require), Profiles will request small amounts of resources for rules that do not require them. This is especially important for rules that are not CPU intensive, as they take far less time to run. This means we can run more jobs at once, and our jobs will finish faster.

### 1. Using the Default Profile
If you are using SLURM, a default profile was downloaded with this repository (under `cluster`). To use it, simply perform the following:
```bash
snakemake --profile cluster
```

This will use the information provided within the `cluster/config.yaml` file to submit jobs to the scheduler.


### 2. Modifying the Profile
Several details of the profile can be modified. Details of the profile are as follows:
- `slurm_account=helikarlab`: If you are not a part of our amazing HelikarLab team, you must change this to your account name. Your account name can be identified by executing the following command on the cluster: `sacctmgr show user $USER accounts`. The account name is the second column of output.
- `restart-times: 0`: This defines the number of times Snakemake shoul

#### 3.1. `config.yaml`
The `config.yaml` is located at `~/.config/snakemake/slurm/config.yaml`. This file contains the configuration for the profile.
Open the file on the cluster, delete its contents, and paste the following:
```yaml
# Default Values
restart-times: 3
jobscript: "slurm-jobscript.sh"
cluster: "slurm-submit.py"
cluster-status: "slurm-status.py"
max-status-checks-per-second: 10
local-cores: 1
latency-wait: 60

# User-modified settings
jobs: 100
printshellcmds: True
max-jobs-per-second: 10  # Default is 1 job per second

# DO NOT MODIFY THESE VALUES
# Snakemake will break if these are changed from their current setting
use-conda: True
conda-frontend: mamba
```

|              Component              |                                                                           Description                                                                            |
|:-----------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|         `restart-times: 3`          |                                                            The number of times to retry a failed job                                                             |
|  `jobscript: "slurm-jobscript.sh"`  |                                                                       The jobscript to use                                                                       |
|    `cluster: "slurm-submit.py"`     |                                                               The cluster submission script to use                                                               |
| `cluster-status: "slurm-status.py"` |                                                                 The cluster status script to use                                                                 |
| `max-status-checks-per-second: 10`  |                                                          The maximum number of status checks per second                                                          |
|          `local-cores: 1`           |                                                            The number of cores to use for local jobs                                                             |
|         `latency-wait: 60`          |                                                The number of seconds to wait before checking the status of a job                                                 |
|             `jobs: 100`             |                                                            The maximum number of jobs to run at once                                                             |
|       `printshellcmds: True`        |                                                            Whether or not to print the shell commands                                                            |
|      `max-jobs-per-second: 10`      |                                                         The maximum number of jobs to submit per second                                                          |
|          `use-conda: True`          |        <span style="color: blue">DO NOT MODIFY</span><br>Whether or not to use Conda to manage dependencies.<br>Modifying this value will break Snakemake        |
|       `conda-frontend: mamba`       | <span style="color: blue">DO NOT MODIFY</span><br>Use mamba when Snakemake installs environments.<br>Modifying this value will significantly slow down Snakemake |


#### 3.2. `cluster_config.yaml`
The `cluster_config.yaml` is located at `~/.config/snakemake/slurm/cluster_config.yaml`. This file contains default configuration values for Slurm to use when submitting jobs to the cluster.
Open the file on the cluster, delete its contents, and paste the following:
```yaml
__default__ :
   job-name  : "{rule}.{wildcards}"
   ntasks    : "1"
   cpus-per-task : "{threads}"
   nodes     : "1"
   output : "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
   error  : "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
```

|    Component    |                                     Description                                      |
|:---------------:|:------------------------------------------------------------------------------------:|
|   `job-name`    |                             The default name of the job                              |
|    `ntasks`     |                          The default number of tasks to run                          |
| `cpus-per-task` |       The default number of cores to use.<br>It is set by each Snakemake rule        |
|     `nodes`     | The number of cluster nodes to request.<br>This should most likely never be changed. |
|    `output`     |                         How (and where) to save output logs                          |
|     `error`     |                          How (and where) to save error logs                          |


## Workflow Configuration

When the workflow was first downloaded (in the [download section](#fastq_download.html)), a `snakemake_config.yaml` file was downloaded as well. Open this file and modify the values to your needs.

To make the `BED_FILE`, `RRNA_INTERVAL_LIST`, and `REF_FLAT_FILE`, see slides 3 and 4 in [this Google Slides][https://docs.google.com/presentation/d/1gxlxbIObhxitgrPLp7lByYFwrFhEdvEm4mILmygAATY/edit#slide=id.g111a3589bd3_0_0] presentation, with examples for the Human Reference Genome

The most up-to-date version of this file can be found [here](https://github.com/HelikarLab/FastqToGeneCounts/blob/436b87c7f40278e918c0ea4b42180243bb84b1d7/snakemake_config.yaml).

### `MASTER_CONTROL`
The master contol file (found under `controls/master_control.csv`) is a CSV file consisting of the following columns:
- SRR codes
- Tissue names + study number, replicate number, and (if applicable) a run number
- Library layouts
    - `PE` for Paired-End
    - `SE` for Single-End
- Library preparation column
    - `total` for total RNA seq
    - `mrna` for PolyA/mRNA RNA seq

If you have `PERFORM_PREFETCH` set to `False` in the `snakemake_config.yaml` file, do not need to modify the master control file. This assumes you are providing the `sra` files yourself. An example of this file is as follows:

{% include warning.html content="The header line should NOT be included in your file" %}

|     SRR     | Tissue/Study/Replicate/Run | Library Layout | Library Prep |
|:-----------:|:--------------------------:|:--------------:|:------------:|
| SRR7647658  |        naiveB_S1R1         |       PE       |     mrna     |
| SRR7647700  |        naiveB_S1R2         |       PE       |     mrna     |
| SRR7647769  |       naiveB_S2R1r1        |       PE       |     mrna     |
| SRR7647808  |       naiveB_S2R1r2        |       PE       |     mrna     |
| SRR5110334  |        naiveB_S3R1         |       SE       |    total     |
| SRR5110338  |        naiveB_S3R2         |       SE       |    total     |
| SRR10408536 |       m2Macro_S1R1r1       |       SE       |    total     |
| SRR10408537 |       m2Macro_S1R1r2       |       SE       |    total     |
| SRR10408538 |       m2Macro_S1R1r3       |       SE       |    total     |
| SRR10408539 |        m2Macro_S2R1        |       SE       |     mrna     |
| SRR10408540 |        m2Macro_S2R2        |       SE       |     mrna     |
| SRR10408541 |        m2Macro_S2R3        |       SE       |     mrna     |


### `DUMP_FASTQ_FILES`
This options is only required if you have set `PERFORM_PREFETCH` to `False`. It is the location at which your input `.fatsq.gz` files are located

### `ROOTDIR`
The relative file path where results should be placed. This is most likely going to be under a `/work` folder. The default value is `results`, which will place results in the `results` folder in the current directory.

### `REF_FLAT_FILE`
The path to a `refFlat` file for your reference genome. This can be made using the following format:
```bash
# Download the gtf to refFlat converter
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

# Add execution permissions
chmod =rwx,g+s ./gtfToGenePred

# Execute the gtf to refFlat converter
# The `genome/Homo_sapiens.GRCh38.105.gtf` is the path to your gtf file
# The last argument, `refFlat.tmp.txt` is the output filename
./gtfToGenePRed -genePredExt -geneNameAsName2 genome/Homo_sapiens.GRCh38.105.gtf refFlat.tmp.txt

# Modify values so Picard is able to parse the refFlat file correctly
# `refFlat.tmp.txt` is the output of the previous command
# `genome/refFlat_GRCh38.105.txt` is the path (and file name) you would like to save results to
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > genome/refFlat_GRCh38.105.txt

# Remove the temporary refFlat file
rm refFlat.tmp.txt
```

### `RRNA_INTERVAL_LIST`
THe path to a ribosomal interval list built from the GTF file for Picard's GetRNASeqMetrics command. This finds rRNA transcript quantities.

The [`riboInt.sh`](https://github.com/HelikarLab/FastqToGeneCounts/blob/436b87c7f40278e918c0ea4b42180243bb84b1d7/riboInt.sh) file was downloaded with the pipeline. You should modify the values so they satisfy your folder paths.

### `BED_FILE`
The path to a BED file for RSeQC, also built from the `GTF_FILE`, which corresponds to your reference genome. This can be made using the following:

```bash
# Change directories into your `genome` directory
cd genome

# Download the BED file
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE.v38.bed.gz/download

# Change the name to a readable name
mv download hg38_GENCODE.v38.bed.gz

# Unzip the file
gunzip hg38_GENCODE.v38.bed.gz

# Remove quotation marks from exon positions
sed -i 's/"//g' hg38_GENCODE.v38.bed

# Remove "chr" from chromosome indices
sed -i 's/chr//g' hg38_GENCODE.v38.bed
```

### `PERFORM_TRIM`
Should trimming of reads be performed? `True` or `False`

### `PERFORM_SCREEN`
Screen against genomes of common contaminants?<br>
`True` or `False`<br>
The current contaminants screened against are:
- [Arabidopsis](https://en.wikipedia.org/wiki/Arabidopsis)
- [Drosophila](https://en.wikipedia.org/wiki/Drosophila)
- [*Escherichia coli*](https://en.wikipedia.org/wiki/Escherichia_coli)
- [Lambda Phage](https://en.wikipedia.org/wiki/Lambda_phage)
- [Mitochondria](https://en.wikipedia.org/wiki/Mitochondrion)
- [Mouse](https://en.wikipedia.org/wiki/House_mouse)
- [PhiX Bacteriophage](https://en.wikipedia.org/wiki/Phi_X_174)
- [Brown Rat](https://en.wikipedia.org/wiki/Brown_rat)
- [rRNA](https://en.wikipedia.org/wiki/Ribosomal_RNA) (specifically, GRCm38 rRNA)
- [Vectors](https://en.wikipedia.org/wiki/Vector_(molecular_biology))
- [*Caenorhabditis elegans*](https://en.wikipedia.org/wiki/Caenorhabditis_elegans) (worm)
- [*Saccharomyces cerevisiae*](https://en.wikipedia.org/wiki/Saccharomyces_cerevisiae) (yeast)

### `PERFORM_GET_RNASEQ_METRICS`
Use Picard's getRNASeqMetrics?<br>
`True` or `False`<br>
This required `REF_FLAT_FILE` and `RRNA_INTERVAL_LIST` to be set.

### `PERFORM_PREFETCH`
If you only have SRR codes (from [`MASTER_CONTROL`](#master_control)), this option will download those `sra` files from NCBI.<br>
`True` or `False`<br>

### `PERFORM_GET_INSERT_SIZE`
Get the interval size using Picard?
`True` or `False`

### `GET_FRAGMENT_SIZE`
Get fragment sizes with RSeQC?
`True` or `False`

### `GENOME_SAVE_DIR`
The path to the directory where genome output should be saved to. This should be under your `/work` folder, as large files will be created

### `GENOME_FASTA_FILE`
This is the input genome fasta file that has been previously downloaded.<br>
Most likely located under your `/work` folder.<br>
Can be downloaded using the following:
```bash
# Change directories into your `genome` directory
cd genome

# Download the assembly
# To set the release number, set the following variable
assembly_release=105
wget ftp://ftp.ensembl.org/pub/release-${assembly_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### `GTF_FILE`
This is the input GTF genome file that has also been previously downloaded<br>
Most likely located under your `/work` folder <br>
Can be downloaded the human genome annotation with:
```bash
# Change directories into your `genome` directory
cd genome

# Download the annotations
# To set releases, modify the following variable to the release number
annotation_release=105
wget ftp://ftp.ensembl.org/pub/release-${annotation_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${annotation_release}.gtf.gz
```


You do not need to make any further changes. Snakemake will extract the configuration values you have set up

Once these steps are complete, the workflow should be prepared to execute.<br>
Continue to the next page to execute the workflow
