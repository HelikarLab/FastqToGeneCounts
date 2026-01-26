# Publication Statistics

This folder contains the scripts used to generate the following statistics:

1. Pearson correlation
2. Spearman correlation
3. Kolmogorov-Smirnov

## Download data
The following SRR IDs must be downloaded before execution (all from `GSE47792`), saved, and unpacked with the `sample` prefix

> [!note]
> You can also chain the [`nf-core/fetchngs`](https://nf-co.re/fetchngs) before the [`nf-core/rnaseq`](`https://nf-co.re/rnaseq`) pipeline to remove manual download and unpacking steps.

```csv
srr,sample
SRR896663,seqc_S1R1r1
SRR896665,seqc_S1R1r2
SRR896667,seqc_S1R1r3
SRR896669,seqc_S1R1r4
SRR896671,seqc_S1R1r5
SRR896673,seqc_S1R1r6
SRR896675,seqc_S1R1r7
SRR896677,seqc_S1R1r8
SRR896679,seqc_S2R1r1
SRR896681,seqc_S2R1r2
SRR896683,seqc_S2R1r3
SRR896685,seqc_S2R1r4
SRR896687,seqc_S2R1r5
SRR896689,seqc_S2R1r6
SRR896691,seqc_S2R1r7
SRR896693,seqc_S2R1r8
SRR896664,seqc_S1R2r1
SRR896666,seqc_S1R2r2
SRR896668,seqc_S1R2r3
SRR896670,seqc_S1R2r4
SRR896672,seqc_S1R2r5
SRR896674,seqc_S1R2r6
SRR896676,seqc_S1R2r7
SRR896678,seqc_S1R2r8
SRR896680,seqc_S2R2r1
SRR896682,seqc_S2R2r2
SRR896684,seqc_S2R2r3
SRR896686,seqc_S2R2r4
SRR896688,seqc_S2R2r5
SRR896690,seqc_S2R2r6
SRR896692,seqc_S2R2r7
SRR896694,seqc_S2R2r8
```

## Download Reference Genome

This pipeline was ran with [Ensembl release 115](https://ftp.ensembl.org/pub/release-115/). The following files must be downloaded for alignment:

1. [Human Primary Assembly](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) 
2. [Human GTF](https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz)
3. [Human Transcriptome](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)

The BED file can be generated with:
```bash
gffread Homo_sapiens*.gtf* -T -o- | awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,$10,$0,$7}' > homo_sapiens.bed
```

To get started, first run the `nf-core/rnaseq` pipeline using nextflow:
```bash
# Generate the configuration file
cat > nfcore.csv <<'EOF'
❯ cat controls/application_note_nfcore.csv
sample,fastq_1,fastq_2,strandedness
seqc_S1R1r1,data/seqc_S1R1r1_1.fastq.gz,data/seqc_S1R1r1_2.fastq.gz,auto
seqc_S1R1r2,data/seqc_S1R1r2_1.fastq.gz,data/seqc_S1R1r2_2.fastq.gz,auto
seqc_S1R1r3,data/seqc_S1R1r3_1.fastq.gz,data/seqc_S1R1r3_2.fastq.gz,auto
seqc_S1R1r4,data/seqc_S1R1r4_1.fastq.gz,data/seqc_S1R1r4_2.fastq.gz,auto
seqc_S1R1r5,data/seqc_S1R1r5_1.fastq.gz,data/seqc_S1R1r5_2.fastq.gz,auto
seqc_S1R1r6,data/seqc_S1R1r6_1.fastq.gz,data/seqc_S1R1r6_2.fastq.gz,auto
seqc_S1R1r7,data/seqc_S1R1r7_1.fastq.gz,data/seqc_S1R1r7_2.fastq.gz,auto
seqc_S1R1r8,data/seqc_S1R1r8_1.fastq.gz,data/seqc_S1R1r8_2.fastq.gz,auto
seqc_S1R2r1,data/seqc_S1R2r1_1.fastq.gz,data/seqc_S1R2r1_2.fastq.gz,auto
seqc_S1R2r2,data/seqc_S1R2r2_1.fastq.gz,data/seqc_S1R2r2_2.fastq.gz,auto
seqc_S1R2r3,data/seqc_S1R2r3_1.fastq.gz,data/seqc_S1R2r3_2.fastq.gz,auto
seqc_S1R2r4,data/seqc_S1R2r4_1.fastq.gz,data/seqc_S1R2r4_2.fastq.gz,auto
seqc_S1R2r5,data/seqc_S1R2r5_1.fastq.gz,data/seqc_S1R2r5_2.fastq.gz,auto
seqc_S1R2r6,data/seqc_S1R2r6_1.fastq.gz,data/seqc_S1R2r6_2.fastq.gz,auto
seqc_S1R2r7,data/seqc_S1R2r7_1.fastq.gz,data/seqc_S1R2r7_2.fastq.gz,auto
seqc_S1R2r8,data/seqc_S1R2r8_1.fastq.gz,data/seqc_S1R2r8_2.fastq.gz,auto
seqc_S1R3r1,data/seqc_S1R3r1_1.fastq.gz,data/seqc_S1R3r1_2.fastq.gz,auto
seqc_S1R3r2,data/seqc_S1R3r2_1.fastq.gz,data/seqc_S1R3r2_2.fastq.gz,auto
seqc_S1R3r3,data/seqc_S1R3r3_1.fastq.gz,data/seqc_S1R3r3_2.fastq.gz,auto
seqc_S1R3r4,data/seqc_S1R3r4_1.fastq.gz,data/seqc_S1R3r4_2.fastq.gz,auto
seqc_S1R3r5,data/seqc_S1R3r5_1.fastq.gz,data/seqc_S1R3r5_2.fastq.gz,auto
seqc_S1R3r6,data/seqc_S1R3r6_1.fastq.gz,data/seqc_S1R3r6_2.fastq.gz,auto
seqc_S1R3r7,data/seqc_S1R3r7_1.fastq.gz,data/seqc_S1R3r7_2.fastq.gz,auto
seqc_S1R3r8,data/seqc_S1R3r8_1.fastq.gz,data/seqc_S1R3r8_2.fastq.gz,auto
seqc_S1R4r1,data/seqc_S1R4r1_1.fastq.gz,data/seqc_S1R4r1_2.fastq.gz,auto
seqc_S1R4r2,data/seqc_S1R4r2_1.fastq.gz,data/seqc_S1R4r2_2.fastq.gz,auto
seqc_S1R4r3,data/seqc_S1R4r3_1.fastq.gz,data/seqc_S1R4r3_2.fastq.gz,auto
seqc_S1R4r4,data/seqc_S1R4r4_1.fastq.gz,data/seqc_S1R4r4_2.fastq.gz,auto
seqc_S1R4r5,data/seqc_S1R4r5_1.fastq.gz,data/seqc_S1R4r5_2.fastq.gz,auto
seqc_S1R4r6,data/seqc_S1R4r6_1.fastq.gz,data/seqc_S1R4r6_2.fastq.gz,auto
seqc_S1R4r7,data/seqc_S1R4r7_1.fastq.gz,data/seqc_S1R4r7_2.fastq.gz,auto
seqc_S1R4r8,data/seqc_S1R4r8_1.fastq.gz,data/seqc_S1R4r8_2.fastq.gz,auto
seqc_S2R1r1,data/seqc_S2R1r1_1.fastq.gz,data/seqc_S2R1r1_2.fastq.gz,auto
seqc_S2R1r2,data/seqc_S2R1r2_1.fastq.gz,data/seqc_S2R1r2_2.fastq.gz,auto
seqc_S2R1r3,data/seqc_S2R1r3_1.fastq.gz,data/seqc_S2R1r3_2.fastq.gz,auto
seqc_S2R1r4,data/seqc_S2R1r4_1.fastq.gz,data/seqc_S2R1r4_2.fastq.gz,auto
seqc_S2R1r5,data/seqc_S2R1r5_1.fastq.gz,data/seqc_S2R1r5_2.fastq.gz,auto
seqc_S2R1r6,data/seqc_S2R1r6_1.fastq.gz,data/seqc_S2R1r6_2.fastq.gz,auto
seqc_S2R1r7,data/seqc_S2R1r7_1.fastq.gz,data/seqc_S2R1r7_2.fastq.gz,auto
seqc_S2R1r8,data/seqc_S2R1r8_1.fastq.gz,data/seqc_S2R1r8_2.fastq.gz,auto
seqc_S2R2r1,data/seqc_S2R2r1_1.fastq.gz,data/seqc_S2R2r1_2.fastq.gz,auto
seqc_S2R2r2,data/seqc_S2R2r2_1.fastq.gz,data/seqc_S2R2r2_2.fastq.gz,auto
seqc_S2R2r3,data/seqc_S2R2r3_1.fastq.gz,data/seqc_S2R2r3_2.fastq.gz,auto
seqc_S2R2r4,data/seqc_S2R2r4_1.fastq.gz,data/seqc_S2R2r4_2.fastq.gz,auto
seqc_S2R2r5,data/seqc_S2R2r5_1.fastq.gz,data/seqc_S2R2r5_2.fastq.gz,auto
seqc_S2R2r6,data/seqc_S2R2r6_1.fastq.gz,data/seqc_S2R2r6_2.fastq.gz,auto
seqc_S2R2r7,data/seqc_S2R2r7_1.fastq.gz,data/seqc_S2R2r7_2.fastq.gz,auto
seqc_S2R2r8,data/seqc_S2R2r8_1.fastq.gz,data/seqc_S2R2r8_2.fastq.gz,auto
EOF

# Generate the params file
cat > nfcore_params.json <<'EOF'
{
  "input": "nfcore.csv",
  "outdir": "./results"
}
EOF

nextflow run nf-core/rnaseq -revision 3.22.2 \
  -resume \
  -profile conda \
  -params-file nfcore_params.json \
  --fasta ./genome/homo_sapiens_release-115/homo_sapiens_release-115_primary_assembly.fa \
  --gtf ./genome/homo_sapiens_release-115/homo_sapiens_release-115.gtf \
  --gene_bed ./genome/homo_sapiens_release-115/homo_sapiens.bed \
  --transcript_fasta ./genome/homo_sapiens_release-115/transcriptome.fa \
  --star_index ./genome/homo_sapiens_release-115/star/ \
  --salmon_index ./genome/homo_sapiens_release-115/salmon/ \
  --max_cpu 22 \
  --skip_bigwig \
  --skip_stringtie \
  --skip_dupradar \
  --skip_qualimap \
  --skip_deseq2_qc \
  --skip_pseudo_alignment \
  --skip_linting
```

A variety of nf-core steps were skipped (`--skip_*`) as they were not needed to validate the output of this pipeline.

## Statistical Tests
### Setup
The `quant.sf` files from the nfcore pipeline and FastqToGeneCounts should be placed into the respective folders `salmon_quant/nfcore` and `salmon_quant/ftgc`

Next, run the associated Python file
```bash
python3 statistics.py
```
