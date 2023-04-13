#!/usr/bin/env bash
# Initially created by Kamil Slowikowski (December 12, 2014)
# Modified by: Arindam Ghosh (July 24, 2019)
# Modified by: Josh Loecker (April 13, 2023)

# Change this variable to the location you would like to save the genome data!
genome_dir="${PWD}/genome"

# Define final output files
primary_assembly_fa="$genome_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
primary_assembly_fai="$genome_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
genome_sizes="$genome_dir/sizes.genome"
genes="$genome_dir/Homo_sapiens.GRCh38.105.gtf"
rRNA_interval_list="$genome_dir/GRCh38.p5.rRNA.interval_list"

# If primary_assembly_fa doesn't exist, show error and quit
if [ ! -e "$primary_assembly_fa" ]; then
  echo "The primary assembly file could not be found\nSearching for:'$primary_assembly_fa'"
fi

# 1. Prepare chromosome sizes file  from fasta sequence if needed.
module load samtools
samtools faidx "$primary_assembly_fa" -o "$primary_assembly_fai"
cut -f1,2 "$primary_assembly_fai" > "$genome_sizes"

# Make an interval_list file suitable for CollectRnaSeqMetrics.jar.
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' "$genome_sizes" | \
    grep -v _ \
>> "$rRNA_interval_list"

# Intervals for rRNA transcripts.
grep 'gene_biotype "rRNA"' $genes | \
    awk '$3 == "gene"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /gene_id "([^"]+)"/ or die "no gene_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> "$rRNA_interval_list"
