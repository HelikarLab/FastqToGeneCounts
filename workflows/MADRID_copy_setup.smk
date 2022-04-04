"""
This file contains rules that are related to copying data to new folders for MADRID
"""

from lib.PerformFunctions import *
import os

rule copy_geneCounts:
    input: rules.star_align.output.gene_table
    output: os.path.join("MADRID_input","{tissue_name}","geneCounts","{sample}","{tissue_name}_{tag}.tab")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        sample=os.path.join("MADRID_input","{tissue_name}","geneCounts","{sample}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
        runtime=1
    shell:
        """
            mkdir -p {params.sample}
            cp {input} {output} || touch {output}
        """

rule copy_strandedness:
    input: rules.get_rnaseq_metrics.output.strand
    output: os.path.join("MADRID_input","{tissue_name}","strandedness","{sample}","{tissue_name}_{tag}_strandedness.txt")

    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        sample=os.path.join("MADRID_input","{tissue_name}","strandedness","{sample}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 50 * attempt,# 0.05 GB
        runtime=1
    shell:
        """
        mkdir -p {params.sample}
        cp {input} {output} || touch {output}
        """

if perform_get_insert_size():
    rule copy_insert_size:
        input: rules.get_insert_size.output.txt
        output: os.path.join("MADRID_input","{tissue_name}","insertSizeMetrics","{sample}","{tissue_name}_{tag}_insert_size.txt")

        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            sample=os.path.join("MADRID_input","{tissue_name}","insertSizeMetrics","{sample}")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
            runtime=1
        shell:
            """
            mkdir -p {params.sample}
            cp {input} {output} || touch {output}
            """

if perform_get_fragment_size():
    rule copy_fragment_size:
        input: rules.get_fragment_size.output
        output: os.path.join("MADRID_input","{tissue_name}","fragmentSizes","{sample}","{tissue_name}_{tag}_fragment_size.txt")

        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            sample=os.path.join("MADRID_input","{tissue_name}","fragmentSizes","{sample}")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 500 * attempt,# 0.5 GB
            runtime=1
        shell:
            """
            mkdir -p {params.sample}
            cp {input} {output} || touch {output}
            """
