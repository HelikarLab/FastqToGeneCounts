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
    print(checkpoint_output)
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

def fastqc_trimmed_reads(wildcards):
    """
    If we are going to trim, return output for rule fastqc_trim
    """
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        return []


rule all:
    input:
        # Generate Genome
        os.path.join(config["ROOTDIR"],config["GENERATE_GENOME"]["GENOME_SAVE_DIR"]),

        # Download SRR codes
        # expand(os.path.join(config["ROOTDIR"], "temp", "rule_complete", "prefetch", "{tissue_name}_{tag}_{srr_code}.complete"), zip, tissue_name=get_tissue_name(), tag=get_tags()(), srr_code=get_srr_data()),

        # dump_fastq
        # This will also request the input of distribute_init_files and prefetch_fastq, without saving their outputs longer than necessary
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),zip,tissue_name=get_tissue_name(),tag=get_tags(),PE_SE=get_PE_SE()),

        # trim reads
        perform_trim,

        # FastQC
        # Untrimed reads (from checkpoint dump_fastq)
        # Trimmed reads (from rule trim)
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),zip,tissue_name=get_tissue_name(),tag=get_tags(),PE_SE=get_PE_SE()),
        fastqc_trimmed_reads,

        # STAR aligner
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab"), zip, tissue_name=get_tissue_name(),tag=get_tags()),

        # MultiQC
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","multiqc","{tissue_name}_multiqc_report.html"), tissue_name=get_tissue_name()),

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
        tissue_name = get_tissue_name()
        tags = get_tags()
        PE_SE = get_PE_SE()
        expanded_files = expand(f"{config['DUMP_FASTQ_FILES']}/{{tissue_name}}_{{tag}}_{{PE_SE}}.fastq.gz", zip, tissue_name=tissue_name, tag=tags, PE_SE=PE_SE)
        for sub_file in expanded_files:
            for sub_tissue, sub_tag, sub_direction in zip(tissue_name, tags, PE_SE):
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

    def get_tag(file_path: str) -> str:
        file_name = os.path.basename(file_path)
        purge_extension = file_name.split(".")[0]
        tag = purge_extension.split("_")[-1]
        return str(tag)
    def get_fastqc_trim_threads(wildcards, input):
        threads = 1
        tag = get_tag(str(input))
        if tag in ["1", "S"]:
            threads = 15
        elif tag == "2":
            threads = 1
        return threads
    def get_fastqc_trim_runtime(wildcards, input, attempt):
        tag = get_tag(str(input))
        runtime = 1
        if tag in ["1", "S"]:
            runtime = 60 * attempt
        elif tag == "2":
            runtime = 5
        return runtime
    rule fastqc_trim:
        input: rules.trim.output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
        params:
            file_two_input=os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_2.fastq.gz"),
            file_two_out=os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_2_fastqc.zip"),
            direction="{PE_SE}"
        threads: get_fastqc_trim_threads
        resources:
            # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
            mem_mb=2048,# 2 GB
            runtime=get_fastqc_trim_runtime
        conda: "envs/fastqc.yaml"
        shell:
            """
            # Process forward reads and reverse reads after trim_galore has finished them
            if [ "{params.direction}" == "1" ]; then
                fastqc {input} --threads {threads} -o $(dirname {output})
                echo "\nFastQC finished $(basename {input}) (1/2)\n"
                fastqc {params.file_two_input} --threads {threads} -o $(dirname {params.file_two_out})
                echo "\nFastQC finished $(basename {params.file_two_input} (2/2)\n"
            elif [ "{params.direction}" == "2" ]; then
                mkdir -p $(dirname {output})
                touch {output}
            elif [ "{params.direction}" == "S" ]; then
                fastqc {input} --threads {threads} -o $(dirname {output})
                echo "\nFastQC finished $(basename {input}) (1/1)\n"
            fi
            """


def collect_star_align_input(wildcards):
    if str(config["PERFORM_TRIM"]).lower() == "true":
        # Have not expanded output from rule trim, need to expand it here
        in_files = expand(rules.trim.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        # already expanding output from dump_fastq, no need to expand it here
        in_files = checkpoints.dump_fastq.get(**wildcards).output

    directions = get_PE_SE()
    grouped_reads = []
    for i, (in_file, direction) in enumerate(zip(in_files, directions)):
        try:
            next_file = in_files[i + 1]
            next_dir = directions[i + 1]
        except:
            if direction == "S":
                grouped_reads.append(in_file)
                continue
            elif direction == "2":
                continue
            else:
                warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
                continue

        if direction == "S":
            grouped_reads.append(in_file)  # "tissue_SXRY_S"
        elif direction == "1" and next_dir == "2":

            next_file = in_files[i + 1]
            if in_file[-9] == next_file[-9]:  # remove _1.fastq.gz to make sure they are same replicate
                both_reads = " ".join([in_file, next_file])  # "tissue_SXRY_1 tissue_SXRY_2"
                grouped_reads.append(both_reads)
            else:
                warnings.warn(f"{in_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {in_file} and one end of {next_file} failed to download")

        elif direction == "1" and not next_dir == "2":
            warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
        elif direction == "2":
            continue
        else:
            warnings.warn(f"{in_file} not handled, unknown reason!")

    """
    We need to return a string, or list of strings. If we return "grouped_reads" directly, some values within are not actually valid files, such as:
        ["results/data/naiveB/naiveB_S1R1_1.fastq.gz results/data/naiveB/naiveB_S1R1_2.fastq.gz", "results/data/naiveB/naiveB_S1R2_S.fastq.gz"]
    Index 0 is taken literally, as a string to a file location. Thus, it does not exist
    Because of this, we are going to filter through each input file and return it if it matches our desired tissue_name and tag
    This is much like what was done in the function get_dump_fastq_output, located above rule all
    """
    for read in grouped_reads:
        if wildcards.tissue_name in read and wildcards.tag in read:
            return read.split(" ")
def get_direction_from_name(file: str):
    file_name = os.path.basename(file)
    purge_extension = file_name.split(".")[0]
    direction = purge_extension.split("_")[-1]
    return direction
def collect_star_align_input_new(wildcards):
    if str(config["PERFORM_TRIM"]).lower() == "true":
        # Have not expanded output from rule trim, need to expand it here
        in_files = sorted(expand(rules.trim.output,zip,tissue_name=get_tissue_name(),tag=get_tags(),PE_SE=get_PE_SE()))
    else:
        # already expanding output from dump_fastq, no need to expand it here
        in_files = sorted(checkpoints.dump_fastq.get(**wildcards).output)

    grouped_reads = []
    for i, in_file in enumerate(in_files):
        direction = get_direction_from_name(in_file)
        try:
            next_file = in_files[i + 1]
            next_direction = get_direction_from_name(next_file)
        except:
            if direction == "S":
                grouped_reads.append(in_file)
                continue
            elif direction == "2": continue
            else:
                warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
                continue

        if direction == "S":
            grouped_reads.append(in_file)
        elif direction == "1" and next_direction == "2":
            if in_file[:-10] == next_file[:-10]:  # remove _1.fastq.gz to make sure they are same replicate
                both_reads = " ".join([in_file, next_file])
                grouped_reads.append(both_reads)
            else:
                warnings.warn(f"{in_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {in_file} and one end of {next_file} failed to download")

        elif direction == "1" and not next_dir == "2":
            warnings.warn(f"{in_file} expects additional paired-end read! Skipping....")
        elif direction == "2":
            continue
        else:
            warnings.warn(f"{in_file} not handled, unknown reason!")

    """
    We need to return a string, or list of strings. If we return "grouped_reads" directly, some values within are not actually valid files, such as:
        ["results/data/naiveB/naiveB_S1R1_1.fastq.gz results/data/naiveB/naiveB_S1R1_2.fastq.gz", "results/data/naiveB/naiveB_S1R2_S.fastq.gz"]
    Index 0 is taken literally, as a string to a file location. Thus, it does not exist
    Because of this, we are going to filter through each input file and return it if it matches our desired tissue_name and tag
    This is much like what was done in the function get_dump_fastq_output, located above rule all
    """
    for read in grouped_reads:
        if wildcards.tissue_name in read and wildcards.tag in read:
            return read.split(" ")


def get_star_align_runtime(wildcards, input, attempt):
    """
    This function will return the length of time required for star_align to complete X number of reads
    Using 40 threads, it takes ~9 minutes per input file
    Round this value to 20 minutes (in case using fewer threads)
    We are also going to multiply by the attempt that the workflow is on.
    If on the second/third/etc. attempt, double/triple/etc. time is requested
    Return an integer of: len(input) * 20 minutes = total runtime
    """
    # Max time is 7 days (10,080 minutes). Do not let this function return more than this time
    return min(len(input.reads) * 20 * attempt, 10079)


rule star_align:
    input:
        reads=collect_star_align_input_new,
        genome_dir=rules.generate_genome.output.genome_dir,
        genome_file=rules.generate_genome.output.genome_file,
        rule_complete=os.path.join(config["ROOTDIR"],"temp","rule_complete","generate_genome.complete")
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        star_output=os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_ReadsPerGene.out.tab")
    conda: "envs/star.yaml"
    threads: 50
    resources:
        mem_mb=51200,# 50 GB
        runtime=get_star_align_runtime
    shell:
        """
        PREFIX="{config[ROOTDIR]}/data/{params.tissue_name}/aligned_reads/{params.tag}/{params.tissue_name}_{params.tag}_"
        echo "prefix is $PREFIX"
        STAR --runThreadN {threads} \
		--readFilesCommand {config[ALIGN_READS][READ_COMMAND]} \
		--readFilesIn {input.reads} \
		--genomeDir {input.genome_dir} \
		--outFileNamePrefix $PREFIX \
		--outSAMtype {config[ALIGN_READS][OUT_SAM_TYPE]} \
		--outSAMunmapped {config[ALIGN_READS][OUT_SAM_UNMAPPED]} \
		--outSAMattributes {config[ALIGN_READS][OUT_SAM_ATTRIBUTES]} \
		--quantMode {config[ALIGN_READS][QUANT_MODE]}

		mv {params.star_output} {output}
        """

def get_fastqc_output(wildcards):
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return expand(rules.fastqc_trim.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        return expand(rules.fastqc_dump_fastq.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
rule multiqc:
    input:
        get_fastqc_output,
        expand(rules.star_align.output, zip, tissue_name=get_tissue_name(), tag=get_tags()),
        lambda wildcards: expand(os.path.join(config["ROOTDIR"], "data", wildcards.tissue_name, "raw", f"{wildcards.tissue_name}_{{tag}}_{{PE_SE}}.fastq.gz"), zip, tag=get_tags(), PE_SE=get_PE_SE())
    output: os.path.join(config["ROOTDIR"],"data", "{tissue_name}","multiqc","{tissue_name}_multiqc_report.html")
    params:
        # lambda not needed as we have tissue_name as wildcard in output
        tissue_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}"),
        tissue_name="{tissue_name}"
    shell:
        """
        multiqc {params.tissue_directory} --filename {params.tissue_name}_multiqc_report.html --outdir {params.tissue_directory}/multiqc/
        """
