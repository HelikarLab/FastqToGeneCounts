#!/usr/bin/env bash
# make_rRNA.sh
# Kamil Slowikowski
# December 12, 2014
#
# Modified: Arindam Ghosh (July 24, 2019 )
#
#
# Referenc Genome: GRCh38.p5 Ensembl release 84
#
#
# 1. Prepare chromosome sizes file  from fasta sequence if needed.
#
#     ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
module load samtools

samtools faidx genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
cut -f1,2 genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > genome/sizes.genome
#
# 2. Make an interval_list file suitable for CollectRnaSeqMetrics.jar.
#
# Ensembl genes:
#
#   ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
#
#
#
# Picard Tools CollectRnaSeqMetrics.jar:
#
#   https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics



chrom_sizes=genome/sizes.genome

# rRNA interval_list file -------------------------------------------------

# Genes from Ensembl.

genes=genome/Homo_sapiens.GRCh38.105.gtf


# Output file suitable for Picard CollectRnaSeqMetrics.jar.

rRNA=genome/GRCh38.p5.rRNA.interval_list

# Sequence names and lengths. (Must be tab-delimited.)
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' $chrom_sizes | \
    grep -v _ \
>> $rRNA

# Intervals for rRNA transcripts.
grep 'gene_biotype "rRNA"' $genes | \
    awk '$3 == "gene"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /gene_id "([^"]+)"/ or die "no gene_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA