# Benchmarking
This is a directory containing scripts to calculate various benchmark data from the pipeline


## Setup
To run the script, simply download the cell-name `.benchmark` files (from the workflow) into this directory.

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


## Running

Simply run `python3 main.py` to start creating graphs for your benchmarks, which will be placed in the `graphs` folder
