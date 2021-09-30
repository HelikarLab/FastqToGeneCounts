import glob
import os
import subprocess
import re
import csv
import numpy as np
import warnings
import math
configfile: "snakemake_config.yaml"

def get_from_master_config(attribute: str) -> list[str]:
    valid_inputs = ["SRR", "tissue", "tag", "PE_SE"]
    sub_list = ["tissue", "tag"]
    if attribute not in valid_inputs:
        sys.exit(f"\nInvalid attribute input. '{attribute}' is not one of: {valid_inputs}\n")
    else:
        collect_attributes = []
        index_value = valid_inputs.index(attribute)

        # We have to subtract one because "tissue" and "tag" are in the same index, thus the index value in valid_inputs is increased by one
        if index_value >= 2: index_value -= 1

        with open(config["MASTER_CONTROL"], "r") as i_stream:
            reader = csv.reader(i_stream)
            for line in reader:
                # Get the column from master_control we are interested in
                column_value = line[index_value]

                # test if we are looking for "tissue" or "tag", as these two values are located at master_control index 1
                if attribute in sub_list:
                    sub_index = sub_list.index(attribute)
                    split_list = str(line[index_value]).split("_")
                    target_attribute = split_list[sub_index]
                    collect_attributes.append(target_attribute)
                else:
                    collect_attributes.append(line[index_value])
        return collect_attributes


def get_srr_code() -> list[str]:
    """
    Only should be getting SRR values if we are performing prefetch
    """
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return get_from_master_config("SRR")

def get_tissue_name() -> list[str]:
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return get_from_master_config("tissue")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tissue_name

def get_tags() -> list[str]:
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return get_from_master_config("tag")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"],"{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tag

def get_PE_SE() -> list[str]:
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return get_from_master_config("PE_SE")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"],"{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.PE_SE

def get_dump_fastq_output(wildcards):
    """
    rule dump_fastq's output has expand(). This means calling rules.dump_fastq.output will retrieve ALL output, not just a single file
    This function will be used to get a single file as input for rule trim, based on rule trim's expected output file
    Example:
        dump_fastq output: results/data/naiveB/raw/naiveB_S1R1_1.fastq.gz, results/data/naiveB/raw/naiveB_S1R1_2.fastq.gz, results/data/naiveB/raw/naiveB_S1R2_1.fastq.gz, results/data/naiveB/raw/naiveB_S1R2_2.fastq.gz
        trim's requested input: results/data/naiveB/raw/naiveB_S1R2_2.fastq.gz
        we will return: results/data/naiveB/raw/naiveB_S1R2_2.fastq.gz
    :param wildcards:
    :return:
    """
    checkpoint_output = checkpoints.dump_fastq.get(**wildcards).output
    for output in checkpoint_output:
        # Only match tissue_name, tag, and PE_SE
        if wildcards.tissue_name in output and \
                wildcards.tag in output and \
                f"_{wildcards.PE_SE}" in output:
            # Getting directory name from dump_fastq
            directory_name = os.path.dirname(output)

            # Generate our working file name, and working file path
            working_file = f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq.gz"
            working_file_path = os.path.join(directory_name, working_file)

            return working_file_path

def perform_trim(wildcards):
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())

    else:
        return []

rule all:
    input:
        # Generate genome
        os.path.join(config["ROOTDIR"],config["GENERATE_GENOME"]["GENOME_SAVE_DIR"]),

        # Dump Fastq
        expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE()),

        # FastQC Dump Fastq
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE()),

        perform_trim,


rule generate_genome:
    input:
        genome_fasta_file=config["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file=config["GENERATE_GENOME"]["GTF_FILE"]
    output:
        genome_dir=directory(os.path.join(config["ROOTDIR"],config["GENERATE_GENOME"]["GENOME_SAVE_DIR"])),
        genome_file=os.path.join(config["ROOTDIR"],config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],"Genome"),
        rule_complete=touch(os.path.join(config["ROOTDIR"],"temp","rule_complete","generate_genome.complete"))
    params:
        log_file=os.path.join(config["ROOTDIR"],"genome","star","Log.out")
    threads: 50
    resources:
        mem_mb=50000, # 50 GB
        runtime=lambda wildcards, attempt: 90 * attempt
    conda: "envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang {config[GENERATE_GENOME][OVERHANG]}

        mv Log.out {params.log_file}
        """

if str(config["PERFORM_PREFETCH"]).lower() == "true":
    rule distribute_init_files:
        input: config["MASTER_CONTROL"]
        output: os.path.join(config["ROOTDIR"],"controls","init_files","{tissue_name}_{tag}.csv")
        params: id="{tissue_name}_{tag}"
        threads: 1
        resources:
            mem_mb=1536,
            runtime=5
        run:
            # Get lines in master control file
            # Open output for writing
            lines = open(str(input),"r").readlines()
            wfile = open(str(output),"w")
            for line in lines:

                # Only write line if the output file has the current tissue-name_tag (naiveB_S1R1) in the file name
                if params.id in line:
                    wfile.write(line)
            wfile.close()

    rule prefetch:
        input: rules.distribute_init_files.output
        output: data=os.path.join(config["ROOTDIR"], "temp", "prefetch", "{tissue_name}_{tag}", "{srr_code}.sra")
        conda: "../../../../../PycharmProjects/FastqToGeneCounts/envs/SRAtools.yaml"
        threads: 1
        resources:
            mem_mb=10240,
            runtime=30
        shell:
            """
            IFS=","
            while read srr name endtype; do
                # prefetch has a default max size of 20G. Effectively remove this size by allowing downloads up to 1TB to be downloaded
                prefetch $srr --max-size 1024000000000 --output-file {output.data}
            done < {input}
            """

def dump_fastq_input(wildcards):
    """
    Return appropriate input for dump_fastq depending on the state of PERFORM_PREFETCH
    """
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return rules.prefetch.output.data
    else:
        tissue_name, tag, PE_SE = glob_wildcards(f"{config['DUMP_FASTQ_FILES']}/{{tissue_name}}_{{tag}}_{{PE_SE}}.fastq.gz")
        expanded_files = expand(f"{config['DUMP_FASTQ_FILES']}/{{tissue_name}}_{{tag}}_{{PE_SE}}.fastq.gz", zip, tissue_name=tissue_name, tag=tag, PE_SE=PE_SE)
        for sub_file in expanded_files:
            for sub_tissue, sub_tag, sub_direction in zip(tissue_name, tag, PE_SE):
                if f"{sub_tissue}" in sub_file and f"_{sub_tag}" in sub_file:
                    return sub_file
checkpoint dump_fastq:
    input: dump_fastq_input
    output: data = os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz")
    script: "scripts/parallel-fastq-dump.py"

rule fastqc_dump_fastq:
    input: lambda wildcards: checkpoints.dump_fastq.get(**wildcards).output
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params: fastqc_output_name=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    shell:
        """
        mkdir ($dirname {output})
        fastqc {input} -o $(dirname {output})
        mv {params.fastqc_output_name} {output}
        echo "\nFastQC finished for {input} (1/1)\n"
        """

if str(config["PERFORM_TRIM"]).lower() == "true":
    def get_trim_threads(wildcards):
        """
        Trim galore uses 9 threads on single-ended data, and 15 cores for paired-end data.
        """
        threads = 1
        if str(wildcards.PE_SE) == "1": threads = 15
        elif str(wildcards.PE_SE) == "2": threads = 1
        elif str(wildcards.PE_SE) == "S": threads = 4
        return threads
    def get_trim_runtime(wildcards, attempt):
        """
        Trim galore takes more time on single-ended data and the forward read of paired end data
        It only touches the output file of the reverse-read paired-end data
        """
        runtime = 5  # minutes
        if str(wildcards.PE_SE) == "1": runtime = 120 * attempt
        elif str(wildcards.PE_SE) == "2": runtime = 5
        elif str(wildcards.PE_SE) == "S": runtime = 120 * attempt
        return runtime
    rule trim:
        input: get_dump_fastq_output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            output_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads"),
            tissue_name="{tissue_name}",
            tag="{tag}",
            direction="{PE_SE}"
        threads: get_trim_threads
        shell:
            """
            # Only process on forward reads
            if [ "{params.direction}" == "1" ]; then
                trim_galore --paired --cores {threads} -o "{params.output_directory}" "{config[ROOTDIR]}/data/{params.tissue_name}/raw/{params.tissue_name}_{params.tag}_1.fastq.gz" "{config[ROOTDIR]}/data/{params.tissue_name}/raw/{params.tissue_name}_{params.tag}_2.fastq.gz"
                filename1="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/{params.tissue_name}_{params.tag}_1_val_1.fq.gz" # final output paired end trimming forward  
                filename2="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/{params.tissue_name}_{params.tag}_2_val_2.fq.gz" # and reverse
                filerename1="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/trimmed_{params.tissue_name}_{params.tag}_1.fastq.gz" # rename to same with trimmed_ prefix
                filerename2="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/trimmed_{params.tissue_name}_{params.tag}_2.fastq.gz" # again       
                mv $filename1 $filerename1 
                mv $filename2 $filerename2
            # Skip over reverse-reads. Create the output file so snakemake does not complain about the rule not generating output
            elif [ "{params.direction}" == "2" ]; then
                touch {output}
            # Work on single-end reads
            elif [ "{params.direction}" == "S" ]; then
                trim_galore --cores {threads} -o "{params.output_directory}" "{input}"
                filename="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/{params.tissue_name}_{params.tag}_S_trimmed.fq.gz" # final out single end trimming
                filerename="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/trimmed_{params.tissue_name}_{params.tag}_{params.direction}.fastq.gz" # rename, same convention as PE
                mv $filename $filerename
            fi
            """

