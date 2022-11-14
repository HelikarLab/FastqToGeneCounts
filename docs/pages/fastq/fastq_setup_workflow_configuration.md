---
title: Workflow Configuration
sidebar: sidebar
permalink: fastq_setup_workflow_configuration.html
summary: This details the workflow configuration
last_updated: Oct 11, 2022
---

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


You do not need to make any further changes. SnakeMake will extract the configuration values you have set up

Once these steps are complete, the workflow should be prepared to execute.<br>
Continue to the next page to execute the workflow
