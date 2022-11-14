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
mkdir genome
cd genome

# Download the assembly
# To set the release number, set the following variable
assembly_release=105
wget ftp://ftp.ensembl.org/pub/release-$(assembly_release)/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

## GTF File
```bash
# Change directories into your `genome` directory
cd genome

# Download the annotations
# To set releases, modify the following variable to the release number
annotation_release=105
wget ftp://ftp.ensembl.org/pub/release-$(annotation_release)/gtf/homo_sapiens/Homo_sapiens.GRCh38.$(annotation_release).gtf.gz
gunzip Homo_sapiens.GRCh38.${annotation_release}.gtf.gz
```

## Ref Flat File
```bash
# Download the gtf to refFlat converter
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

# Add execution permissions
chmod +x ./gtfToGenePred

# Execute the gtf to refFlat converter
# The `genome/Homo_sapiens.GRCh38.105.gtf` is the path to your gtf file
# The last argument, `refFlat.tmp.txt` is the output filename
./gtfToGenePred -genePredExt -geneNameAsName2 genome/Homo_sapiens.GRCh38.105.gtf refFlat.tmp.txt

# Modify values so Picard is able to parse the refFlat file correctly
# `refFlat.tmp.txt` is the output of the previous command
# `genome/refFlat_GRCh38.105.txt` is the path (and file name) you would like to save results to
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > genome/refFlat_GRCh38.105.txt

# Remove the temporary refFlat file
rm refFlat.tmp.txt
```

## BED File
```bash
# Change directories into your `genome` directory
cd genome

# Download the assembly
# To set the release number, set the following variable
assembly_release=105
wget ftp://ftp.ensembl.org/pub/release-$(assembly_release)/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```


## rRNA Interval List
The path to a ribosomal interval list built from the GTF file for Picard’s GetRNASeqMetrics command. This finds rRNA transcript quantities.

The [`riboInt.sh`](https://github.com/HelikarLab/FastqToGeneCounts/blob/436b87c7f40278e918c0ea4b42180243bb84b1d7/riboInt.sh) file was downloaded with the pipeline.

In theory, this file should fit the current paths set up according to these instructions. However, you should double check the values within satisfy your own setup.
