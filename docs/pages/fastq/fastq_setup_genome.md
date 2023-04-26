---
title: Setting up the Genome Folder
sidebar: sidebar
permalink: fastq_setup_genome.html
summary: This is an overview of how to set up the `genome` folder
last_updated: Oct 11, 2022
---

## Overview
This will overview how to set up the `genome` folder.

## Genome FASTA File
```bash
# Change directories into your `genome` directory
cd FastqToGeneCounts

mkdir genome
cd genome

# Download the assembly
# To set the release number, set the following variable
assembly_release=105
wget ftp://ftp.ensembl.org/pub/release-${assembly_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.*.dna.primary_assembly.fa.gz
```

## GTF File
```bash
# Execute this in the `genome` directory!

# Download the annotations
# To set releases, modify the following variable to the release number
annotation_release=105
wget ftp://ftp.ensembl.org/pub/release-${annotation_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${annotation_release}.gtf.gz
gunzip Homo_sapiens.GRCh38.${annotation_release}.gtf.gz
```

## Ref Flat File
```bash
# Execute this in the `genome` directory!

# Download the gtf to refFlat converter
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

# Add execution permissions
chmod +x ./gtfToGenePred

# Execute the gtf to refFlat converter
# The last argument, `refFlat.tmp.txt` is the output filename
./gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.105.gtf refFlat.tmp.txt

# Modify values so Picard is able to parse the refFlat file correctly
# `refFlat.tmp.txt` is the output of the previous command
# `refFlat_GRCh38.105.txt` is the final output filename
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat_GRCh38.105.txt

# Remove the temporary refFlat file
rm refFlat.tmp.txt
```

## BED File
```bash
# Execute this in the `genome` directory!

# Download the BED file, then set the file name
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE.v38.bed.gz/download
mv download hg38_GENCODE.v38.bed.gz

# Unzip the file
gunzip hg38_GENCODE.v38.bed.gz

# remove (1) quotation marks from exon positions and (2) “chr” from chromosome indices
sed -i 's/"//g' hg38_GENCODE.v38.bed
sed -i 's/chr//g' hg38_GENCODE.v38.bed
```

## rRNA Interval List
The path to a ribosomal interval list built from the GTF file for Picard’s GetRNASeqMetrics command. This finds rRNA transcript quantities.

The [`riboInt.sh`](https://github.com/HelikarLab/FastqToGeneCounts/blob/436b87c7f40278e918c0ea4b42180243bb84b1d7/riboInt.sh) file was downloaded with the pipeline.

In theory, this file should fit the current paths set up according to these instructions. However, you should double check the values within satisfy your own setup.

To run the riboInt.sh file, perform the following:

```bash
# Change directories to the location the pipeline was downloaded; for example:
cd /work/helikarlab/joshl/FastqToGeneCounts
sh riboInt.sh
```
