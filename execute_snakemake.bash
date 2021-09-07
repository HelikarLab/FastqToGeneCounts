#!/bin/bash

#SBATCH --job-name="snakemake"
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 06:00:00
#SBATCH -o "slurm/snakemake.out"
#SBATCH -e "slurm/snakemake.err"

source activate /home/helikarlab/joshl/work/fastq_dump
module load parallel-fastq-dump
module load SRAtoolkit

cd /home/helikarlab/joshl/work/b_cell

date
echo "Starting"
snakemake --cores all
