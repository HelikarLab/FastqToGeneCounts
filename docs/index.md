---
title: Getting Started
sidebar: sidebar
tags: [getting_started]
permalink: index.html
summary: This is an overview of how to get started with FastqToGeneCounts
last_updated: Sept 20, 2022
---

# Welcome!

This is a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that aims to do several things, using as much parallelization as possible.

Given a CSV file containing: SRR codes, target output names, library layouts, and library preparation methods:

1. Generate genome files using [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
2. Download each SRR code in parallel using [prefetch](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)
3. Unpack the `.sra` files using [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump), generating `.fastq.gz` files
4. Optionally trim the resulting `.fastq.gz` files (using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
5. Perform [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the parallel-fastq-dump files, and optionally on the resulting trimmed files
6. Perform STAR align on files from parallel-fastq-dump (or trim) files to the generated genome files
7. Optionally get RNAseqMetrics using [Picard](https://broadinstitute.github.io/picard/)
8. Optionally get insert sizes using Picard
9. Optionally get fragment sizes using [RSeQC](http://rseqc.sourceforge.net/)
10. Perform [MultiQC](https://multiqc.info), using the files from parallel-fastq-dump, FastQC, and STAR aligner
11. Organize a MADRID_inputs file that can be directly interfaced with [our MADRID package](https://github.com/HelikarLab/MADRID) to aid with metabolic drug discovery and repurposing.

This pipeline is primarily designed to interface the [GEO Database](https://www.ncbi.nlm.nih.gov/geo/) with MADRID, and should be run in a high-performance computing cluster, as the memory requirement is quite high to use STAR (about 40GB for the human genome). Even if you do not plan to use MADRID, if your goal is to align fastq files from bulk RNA-seq, perform essential quality control, and output gene counts files from STAR for transcription-based model construction such as Differential Gene Expression Analysis, this pipeline could be of service.
