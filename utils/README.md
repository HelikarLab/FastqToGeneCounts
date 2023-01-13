# Utils

This is a directory for various utilities used throughout the FastqToGeneCounts project.  
The most useful file in this folder is the `benchmark.py` file, which will create a series of plots using data from Snakemake's benchmarking capabilities. Read further on how to create benchmark graphs

# Benchmarking
This is a directory containing scripts to calculate various benchmark data from the pipeline

## Setup
Simply execute `python3 benchmark.py` to create benchmark graphs. The `benchmark.py` file will search through the `benchmarks` directory, one folder "above" this one.

The following is an example of `preB` cell type benchmark information

**Example 1**
```
preB
├── contaminant_screen
│   ├── preB_S1R1_1.benchmark
│   ├── preB_S1R1_2.benchmark
│   ├── preB_S1R2_1.benchmark
│   ├── preB_S1R2_2.benchmark
│   ├── preB_S1R3_1.benchmark
│   ├── preB_S1R3_2.benchmark
│   ├── preB_S1R4_1.benchmark
│   └── preB_S1R4_2.benchmark
├── copy_geneCounts
│   ├── S1_preB_S1R1.benchmark
│   ├── S1_preB_S1R2.benchmark
│   ├── S1_preB_S1R3.benchmark
│   └── S1_preB_S1R4.benchmark
├── copy_insert_size
│   ├── S1_preB_S1R1.benchmark
│   ├── S1_preB_S1R2.benchmark
│   ├── S1_preB_S1R3.benchmark
│   └── S1_preB_S1R4.benchmark
├── copy_strandedness
│   ├── S1_preB_S1R1.benchmark
│   ├── S1_preB_S1R2.benchmark
│   └── S1_preB_S1R3.benchmark
├── distribute_init_files
│   ├── preB_S1R1.benchmark
│   ├── preB_S1R2.benchmark
│   ├── preB_S1R3.benchmark
│   └── preB_S1R4.benchmark
├── fasterq_dump
│   ├── preB_S1R1_1.benchmark
│   ├── preB_S1R1_2.benchmark
│   ├── preB_S1R2_1.benchmark
│   ├── preB_S1R2_2.benchmark
.
.
.
```