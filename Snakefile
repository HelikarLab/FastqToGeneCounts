import glob
import os
import subprocess
import re
import csv
import numpy as np
import warnings

configfile: "snakemake_config.yaml"

def get_tissue_name():
    """
    Looking to return the base filename from the controls/init_files
    Example:
    controls/init_files/naiveB_S1R1.csv
    controls/init_files/naiveB_S1R2.csv
    We would return: ["naiveB", "naiveB"]
    :return:
    """
    tissue_data = []

    with open(config["MASTER_CONTROL"],"r") as rfile:
        reader = csv.reader(rfile)
        for i, line in enumerate(reader):
            id = line[1].split("_")[0]  # naiveB_S1R1 -> naiveB
            pe_se = line[2]

            # append tissue name twice of paired end, allows for naming "_1" and "_2"
            if pe_se == "PE":
                tissue_data.append(id)
                tissue_data.append(id)
            elif pe_se == "SE":
                tissue_data.append(id)

    return tissue_data


def get_tag_data():
    """
    Return tag from cell ID
    Example:
        input: naiveB_S1R1
        output: S1R1
    :return:
    """
    tag_data = []
    with open(config["MASTER_CONTROL"],"r") as rfile:
        reader = csv.reader(rfile)
        for i, line in enumerate(reader):
            tag = line[1].split("_")[-1]
            pe_se = line[2]

            # append tag twice for paired end, allows for naming "_1" and "_2"
            if pe_se == "PE":
                tag_data.append(tag)
                tag_data.append(tag)
            elif pe_se == "SE":
                tag_data.append(tag)
    return tag_data


def get_srr_data():
    """
    Get the SRR information from the master init_file
    Example:
        input:
            SRR14231328,naiveB_S1R1,PE
            SRR14231329,naiveB_S1R2,PE
        output: ["SRR14231329", "SRR14231328"]
    :return:
    """
    srr_data = []
    with open(config["MASTER_CONTROL"],"r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            srr = line[0]
            pe_se = line[2]

            if pe_se == "PE":
                srr_data.append(srr)
                srr_data.append(srr)
            elif pe_se == "SE":
                srr_data.append(srr)
    return srr_data


def get_PE_SE_Data():
    """
    This function will read from the config[MASTER_CONTROL] file and return the paired_end or single_end variable
    Example:
        input:
            SRR14231328,naiveB_S1R1,PE
            SRR14231329,naiveB_S1R2,SE
        output:
            ["_1", "_2", "_s"]
            # PE,   PE,   SE
    :return:
    """
    pe_se_data = []
    with open(config["MASTER_CONTROL"],"r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            pe_se = line[2]
            if pe_se == "PE":
                pe_se_data.append("1")
                pe_se_data.append("2")
            elif pe_se == "SE":
                pe_se_data.append("S")
    return pe_se_data


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
            working_file_path = os.path.join(directory_name,working_file)

            return working_file_path


rule all:
    input:
        # Generate Genome
        os.path.join(config["ROOTDIR"],config["GENERATE_GENOME"]["GENOME_SAVE_DIR"]),

        # dump_fastq
        # This will also request the input of distribute_init_files and prefetch_fastq, without saving their outputs longer than necessary
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),zip,tissue_name=get_tissue_name(),tag=get_tag_data(),PE_SE=get_PE_SE_Data()),

        # Trim Reads
        # May not need to include this, as it is dynamic depending on what STAR aligner needs
        # If config["PERFORM_TRIM"] == False, including this as input will cause an error
        # expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tag_data(), PE_SE=get_PE_SE_Data()),

        # FastQC
        # Untrimed reads (i.e. from dump fastq)
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),zip,tissue_name=get_tissue_name(),tag=get_tag_data(),PE_SE=get_PE_SE_Data()),
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip"), zip, tissue_name=get_tissue_name(), tag=get_tag_data(), PE_SE=get_PE_SE_Data()),

        # STAR aligner
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_ReadsPerGene.out.tab"),zip,tissue_name=get_tissue_name(),tag=get_tag_data()),
        #directory(os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads"))

        # MultiQC
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","multiqc","{tissue_name}_multiqc_report.html"),tissue_name=get_tissue_name())

#TODO: convert this input to snakemake's HTTP download
# https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#read-only-web-http-s
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
    threads: workflow.cores * 0.9
    resources:
        mem_mb = 50000, # 50 GB
        runtime = 45    # 45 minutes
    conda: "envs/STAR.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang {config[STAR][GENERATE_GENOME][OVERHANG]}

        mv Log.out {params.log_file}
        """

rule distribute_init_files:
    input: config["MASTER_CONTROL"]
    output: temp(os.path.join(config["ROOTDIR"],"controls","init_files","{tissue_name}_{tag}.csv"))
    params: id="{tissue_name}_{tag}"
    threads: 1
    resources:
        mem_mb = 1536,  # 1.5 GB
        runtime = 5    # 5 minutes
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

rule prefetch_fastq:
    input: rules.distribute_init_files.output
    output: data=temp(os.path.join(config["ROOTDIR"],"temp","prefetch","{tissue_name}_{tag}","{srr_code}","{srr_code}.sra"))
    envmodules: "SRAtoolkit/2.10"
    threads: 1
    resources:
        mem_mb = 10240, # 10 GB
        runtime = 30    # 30 minutes
    shell:
        """
        IFS=","
        while read srr name endtype; do
            # prefetch has a default max size of 20G. Effectively remove this size by allowing downloads up to 1TB to be downloaded
            prefetch $srr --max-size 1024000000000 --output-file {output.data}
        done < {input}
        """


def get_dump_fastq_runtime(wildcards, input, attempt):
    """
    This function will dynamicall return the length of time requested by checkpoint dump_fastq
    Using 40 threads, it takes approximately 5 minutes per input file
    We are going to round this up to 20 minutes (= 1200 seconds) to be extremely safe

    The 'attempt' input is a snakemake global variable.
    If this workflow fails when using the profile option "restart-times" or the command line option "--restart-times"
        the workflow will be restarted X many times. When doing so, we will increase the amount of time requested
        on each subsequent attempt.
        i.e. attempt 2 will double the amount of time requested for this rule, attempt 3 will triple the amount of time requested

    :param wildcards: wildcard input from checkpoint dump_fastq
    :param input: all input files
    :param attempt: the attempt number of the run
    :return: integer, length of input * 20 minutes
    """
    return len(input) * 1200 * attempt

"""
Shifted this section to using "scripts" and "conda" because it allows us to package everything we need within the script itself
We are able to download the parallel-fastq-dump pacakge from Anaconda where-ever we are, and do not depend on a cluster having it installed 
"""
checkpoint dump_fastq:
    input: data = expand(rules.prefetch_fastq.output.data,zip,tissue_name=get_tissue_name(),tag=get_tag_data(),srr_code=get_srr_data())
    output: data = expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tag_data(), PE_SE=get_PE_SE_Data())
    threads: workflow.cores * 0.9  # max threads
    conda: "envs/parallel-fastq-dump.yaml"
    resources:
        mem_mb=20480,  # 20 GB
        runtime=get_dump_fastq_runtime
    script:
        "scripts/parallel-fastq-dump.py"


rule fastqc_dump_fastq:
    input: get_dump_fastq_output
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params: outdir=os.path.dirname(os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","{tissue_name}_{tag}_{PE_SE}_fastqc.zip"))
    threads: 5
    resources:
        # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
        mem_mb = 2048,# 2 GB
        runtime = 60  # 60 minutes
    envmodules: "fastqc"
    shell:
        """
        fastqc {input} --threads {threads} -o {params.outdir}
        """

if str(config["PERFORM_TRIM"]).lower() == "true":
    rule trim:
        input: get_dump_fastq_output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            output_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads"),
            working_file="{tissue_name}_{tag}_{PE_SE}.fastq.gz",
            dump_fastq_output_dir=get_dump_fastq_output,
            tissue_name="{tissue_name}",
            tag="{tag}",
            direction="{PE_SE}"
        envmodules: "trim_galore/0.6"
        threads: 1
        # TODO: Get resources for this rule
        shell:
            """
            # Only process on forward reads
            if [ "{params.direction}" == "1" ]; then
                trim_galore --paired -o "{params.output_directory}" "{config[ROOTDIR]}/data/{params.tissue_name}/raw/{params.tissue_name}_{params.tag}_1.fastq.gz" "{config[ROOTDIR]}/data/{params.tissue_name}/raw/{params.tissue_name}_{params.tag}_2.fastq.gz"
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
                trim_galore -o "{params.output_directory}" "{input}"
                filename="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/{params.tissue_name}_{params.tag}_S_trimmed.fq.gz" # final out single end trimming
                filerename="{config[ROOTDIR]}/data/{params.tissue_name}/trimmed_reads/trimmed_{params.tissue_name}_{params.tag}_{params.direction}.fastq.gz" # rename, same convention as PE
                mv $filename $filerename
            fi

            """

    rule fastqc_trim:
        input: rules.trim.output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
        threads: 5
        resources:
            # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
            mem_mb=2048,# 2 GB
            runtime=60  # 60 minutes
        envmodules: "fastqc"
        shell:
            """
            fastqc {input} --threads {threads} -o {output}
            """

def collect_star_align_input(wildcards):
    if str(config["PERFORM_TRIM"]).lower() == "true":
        # Have not expanded output from rule trim, need to expand it here
        in_files = expand(rules.trim.output,zip,tissue_name=get_tissue_name(),tag=get_tag_data(),PE_SE=get_PE_SE_Data())
    else:
        # already expanding output from dump_fastq, no need to expand it here
        in_files = checkpoints.dump_fastq.get(**wildcards).output

    directions = get_PE_SE_Data()
    grouped_reads = []
    for i, (in_file, direction) in enumerate(zip(in_files,directions)):
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

def get_star_align_runtime(wildcards, input, attempt):
    """
    This function will return the length of time required for star_align to complete X number of reads

    Using 40 threads, it takes ~9 minutes per input file
    Round this value to 20 minutes (in case using fewer threads)

    We are also going to multiply by the attempt that the workflow is on.
    If on the second/third/etc. attempt, double/triple/etc. time is requested

    Return an integer of: len(input) * 1200 seconds = total runtime
    :param wildcards:
    :param input:
    :return:
    """
    return len(input.reads) * 1200 * attempt
rule star_align:
    input:
        reads=collect_star_align_input,
        genome_dir=rules.generate_genome.output.genome_dir,
        genome_file=rules.generate_genome.output.genome_file,
        rule_complete=os.path.join(config["ROOTDIR"],"temp","rule_complete","generate_genome.complete")
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_ReadsPerGene.out.tab")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}"
    envmodules: "star/2.7"
    threads: workflow.cores * 0.90
    resources:
        mem_mb=51200,# 50 GB
        runtime=get_star_align_runtime
    shell:
        """
        PREFIX="{config[ROOTDIR]}/data/{params.tissue_name}/aligned_reads/{params.tag}/{params.tissue_name}_{params.tag}_"
        echo "prefix is $PREFIX"
        STAR --runThreadN {threads} \
		--readFilesCommand {config[STAR][ALIGN_READS][READ_COMMAND]} \
		--readFilesIn {input.reads} \
		--genomeDir {input.genome_dir} \
		--outFileNamePrefix $PREFIX \
		--outSAMtype {config[STAR][ALIGN_READS][OUT_SAM_TYPE]} \
		--outSAMunmapped {config[STAR][ALIGN_READS][OUT_SAM_UNMAPPED]} \
		--outSAMattributes {config[STAR][ALIGN_READS][OUT_SAM_ATTRIBUTES]} \
		--quantMode {config[STAR][ALIGN_READS][QUANT_MODE]}
        """

rule multiqc:
    input:
        # We are using "lambda wildcards" here so we do not have to use a function that contains only "return checkpoints.dump_fastq.get(**wildcards).output"
        lambda wildcards: checkpoints.dump_fastq.get(**wildcards).output,
        expand(rules.fastqc_dump_fastq.output,zip,tissue_name=get_tissue_name(),tag=get_tag_data(),PE_SE=get_PE_SE_Data()),
        expand(rules.star_align.output,zip,tissue_name=get_tissue_name(),tag=get_tag_data())
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","multiqc","{tissue_name}_multiqc_report.html")
    params:
        # lambda not needed as we have tissue_name as wildcard in output
        tissue_directory = os.path.join(config["ROOTDIR"], "data", "{tissue_name}"),
        tissue_name = "{tissue_name}"
    envmodules: "multiqc/py37/1.8"
    threads: 1
    resources:
        # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
        mem_mb=1024,  # 1 GB
        runtime=10    # 10 minutes
    shell:
        """
        multiqc {params.tissue_directory} --filename {params.tissue_name}_multiqc_report.html --outdir {params.tissue_directory}/multiqc/
        """
