import glob
import os
import subprocess
import re
import csv
import numpy as np
import warnings
import math
configfile: "snakemake_config.yaml"

def perform_trim():
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return True
    else:
        return False
def perform_prefetch():
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return True
    else:
        return False

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
        i_stream = open(config["MASTER_CONTROL"], "r")
        reader = csv.reader(i_stream)
        for line in reader:

            # Get the column from master_control we are interested in
            column_value = line[index_value]
            PE_SE_value = line[2]  # get PE or SE

            # test if we are looking for "tissue" or "tag", as these two values are located at master_control index 1
            if attribute in sub_list:
                sub_index = sub_list.index(attribute)
                split_list = str(line[index_value]).split("_")

                # We must append the target attribute twice if it is paired end, once if it is single end
                if PE_SE_value == "PE":
                    target_attribute = [split_list[sub_index], split_list[sub_index]]
                else:  # Single end
                    target_attribute = [split_list[sub_index]]

            elif attribute == "PE_SE":
                # We must append the target attribute twice if it is paired end, once if it is single end
                if column_value == "PE":
                    target_attribute = ["1", "2"]
                else:  # Single end
                    target_attribute = ["S"]

            else:
                if PE_SE_value == "PE":
                    target_attribute = [line[index_value], line[index_value]]
                else:
                    target_attribute = [line[index_value]]

            collect_attributes += target_attribute

        i_stream.close()
        return collect_attributes

def get_srr_code() -> list[str]:
    """
    Only should be getting SRR values if we are performing prefetch
    """
    if perform_trim():
        return get_from_master_config("SRR")

def get_tissue_name() -> list[str]:
    if perform_prefetch():
        return get_from_master_config("tissue")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tissue_name

def get_tags() -> list[str]:
    if perform_prefetch():
        return get_from_master_config("tag")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tag

def get_PE_SE() -> list[str]:
    if perform_prefetch():
        return get_from_master_config("PE_SE")
    else:
        fastq_input = glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
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

def perform_trim_rule(wildcards):
    if perform_trim():
        return expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        return []

def fastqc_trimmed_reads(wildcards):
    """
    If we are going to trim, return output for rule fastqc_trim
    """
    if perform_trim():
        return expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        return []

def perform_dump_fastq(wildcards):
    if perform_prefetch():
        data = expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
        return data
    else:
        return []


rule all:
    input:
        # Generate Genome
        os.path.join(config["ROOTDIR"], config["GENERATE_GENOME"]["GENOME_SAVE_DIR"]),

        # dump_fastq
        perform_dump_fastq,

        # trim reads
        perform_trim_rule,

        # FastQC
        # Untrimed reads (from checkpoint dump_fastq)
        # Trimmed reads (from rule trim)
        expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE()),
        fastqc_trimmed_reads,

        # STAR aligner
        expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.tab"), zip, tissue_name=get_tissue_name(), tag=get_tags()),

        # MultiQC
        expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "multiqc", "{tissue_name}_multiqc_report.html"), tissue_name=get_tissue_name()),

rule generate_genome:
    input:
        genome_fasta_file=config["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file=config["GENERATE_GENOME"]["GTF_FILE"]
    output:
        genome_dir=directory(os.path.join(config["ROOTDIR"], config["GENERATE_GENOME"]["GENOME_SAVE_DIR"])),
        genome_file=os.path.join(config["ROOTDIR"], config["GENERATE_GENOME"]["GENOME_SAVE_DIR"], "Genome"),
        rule_complete=touch(os.path.join(config["ROOTDIR"], "temp", "rule_complete", "generate_genome.complete"))
    params:
        log_file=os.path.join(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"], "Log.out")
    threads: 40
    resources:
        mem_mb=50000, # 50 GB
        runtime=lambda wildcards, attempt: 120 * attempt
    conda: "envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99

        mv Log.out {params.log_file}
        """

def get_dump_fastq_runtime(wildcards, input, attempt):
    """
    This function will dynamicall return the length of time requested by checkpoint dump_fastq
    Using 40 threads, it takes approximately 5 minutes per input file
    We are going to round this up to 10 minutes to be safe
    The 'attempt' input is a snakemake global variable.
    If this workflow fails when using the profile option "restart-times" or the command line option "--restart-times"
        the workflow will be restarted X many times. When doing so, we will increase the amount of time requested
        on each subsequent attempt.
        i.e. attempt 2 will double the amount of time requested for this rule, attempt 3 will triple the amount of time requested
    We are dividing the input by 2 because duplicates are input, one for the forward read and one for the reverse read
    :param wildcards: wildcard input from checkpoint dump_fastq
    :param input: all input files
    :param attempt: the attempt number of the run
    :return: integer, length of input * 20 minutes
    """
    # Max time is 10,080 minutes (7 days), do not let this function return more than that amount of time
    return 30 * attempt
if perform_prefetch():
    rule distribute_init_files:
        input: config["MASTER_CONTROL"]
        output: os.path.join(config["ROOTDIR"], "controls", "init_files", "{tissue_name}_{tag}.csv")
        params: id="{tissue_name}_{tag}"
        threads: 1
        resources:
            mem_mb=1536,
            runtime=5
        run:
            # Get lines in master control file
            # Open output for writing
            lines = open(str(input), "r").readlines()
            wfile = open(str(output), "w")
            for line in lines:

                # Only write line if the output file has the current tissue-name_tag (naiveB_S1R1) in the file name
                if params.id in line:
                    wfile.write(line)
            wfile.close()

    rule prefetch:
        input: rules.distribute_init_files.output
        output: os.path.join(config["ROOTDIR"], "temp", "prefetch", "{tissue_name}_{tag}", "{srr_code}.sra")
        conda: "envs/SRAtools.yaml"
        threads: 1
        resources:
            mem_mb=10240,
            runtime=lambda wildcards, attempt: 30 * attempt
        shell:
            """
            IFS=", "
            while read srr name endtype; do
                # prefetch has a default max size of 20G. Effectively remove this size by allowing files up to 1TB to be downloaded
                prefetch $srr --max-size 1024000000000 --output-file {output}
            done < {input}
            """

    def dump_fastq_input(wildcards):
        output_files = expand(rules.prefetch.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), srr_code=get_srr_code())
        for file in output_files:
            if (wildcards.tissue_name in file) and (wildcards.tag in file):
                return file
    def get_dump_fastq_threads(wildcards):
        """Get threads for dump fastq"""
        threads = 1
        if str(wildcards.PE_SE) in ["1", "S"]: threads = 40
        elif str(wildcards.PE_SE) == "2": threads = 1
        return threads
    def get_dump_fastq_srr_code(wildcards):
        """Get SRR codes corresponding to dump_fastq output"""
        lines = open(config["MASTER_CONTROL"], "r").readlines()
        if wildcards.PE_SE in ["1", "2"]: direction = "PE"
        else: direction = "SE"
        for line in lines:
            if f"{wildcards.tissue_name}_{wildcards.tag},{direction}" in line:
                srr_code = line.split(",")[0]
                return srr_code

    checkpoint dump_fastq:
        input: dump_fastq_input
            # expand(rules.prefetch.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), srr_code=get_srr_code())
        output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            srr_code=get_dump_fastq_srr_code
        threads: get_dump_fastq_threads
        conda: "envs/SRAtools.yaml"
        resources:
            mem_mb=20480,  # 20 GB
            runtime=get_dump_fastq_runtime
        shell:
            """
            output_dir="$(dirname {output})"
            if [ "{wildcards.PE_SE}" == "1" ]; then
                parallel-fastq-dump --sra-id {input} --threads {threads} --outdir "$output_dir" --gzip --split-files
                
                mv "$output_dir/{params.srr_code}_1.fastq.gz" "$output_dir/{wildcards.tissue_name}_{wildcards.tag}_1.fastq.gz"
                mv "$output_dir/{params.srr_code}_2.fastq.gz" "$output_dir/{wildcards.tissue_name}_{wildcards.tag}_2.fastq.gz"
                 
            elif [ "{wildcards.PE_SE}" == "2" ]; then
                touch {output}
            elif [ "{wildcards.PE_SE}" == "S" ]; then
                parallel-fastq-dump --sra-id {input} --threads {threads} --outdir "$output_dir" --gzip
                
                mv "$output_dir/{params.srr_code}.fastq.gz" "$output_dir/{wildcards.tissue_name}_{wildcards.tag}_S.fastq.gz"
            fi
            """


def get_tag(file_path: str) -> str:
    file_name = os.path.basename(file_path)
    purge_extension = file_name.split(".")[0]
    tag = purge_extension.split("_")[-1]
    return str(tag)
def get_fastqc_threads(wildcards, input):
    threads = 1
    tag = get_tag(str(input))
    if tag in ["1", "S"]:
        threads = 15
    elif tag == "2":
        threads = 1
    return threads
def get_fastqc_runtime(wildcards, input, attempt):
    tag = get_tag(str(input))
    runtime = 1
    if tag in ["1", "S"]:
        runtime = 150 * attempt  # 2.5 hours
    elif tag == "2":
        runtime = 5
    return runtime
def fastqc_dump_fastq_input(wildcards):
    if perform_prefetch():
        return checkpoints.dump_fastq.get(**wildcards).output
    else:
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                if (wildcards.tissue_name in file) and (wildcards.tag in file) and (f"_{wildcards.PE_SE}" in file):
                    return os.path.join(path, file)
rule fastqc_dump_fastq:
    input: fastqc_dump_fastq_input
    output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "untrimmed_reads", "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params:
        fastqc_output_name=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "untrimmed_reads", "{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        file_two_input=os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_2.fastq.gz"),
        file_two_out=os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_2_fastqc.zip"),
        direction="{PE_SE}"
    threads: get_fastqc_threads
    conda: "envs/fastqc.yaml"
    resources:
        mem_mb=2048,  # 2 GB
        runtime=get_fastqc_runtime    # 150 minutes * attempt number
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"
        
        if [ "{params.direction}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"
            printf "FastQC finished $(basename {input}) (1/2)\n\n"
            
            fastqc {params.file_two_input} --threads {threads} -o "$output_directory"
            printf "FastQC finished $(basename {params.file_two_input}) (2/2)\n\n"
        elif [ "{params.direction}" == "2" ]; then
            touch {output}
        elif [ "{params.direction}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"
            printf "FastQC finished $(basename {input}) (1/1)\n\n"
        fi
        """

if perform_trim():
    def get_trim_input(wildcards):
        if perform_prefetch():
            return checkpoints.dump_fastq.get(**wildcards).output
        else:
            fastq_gz_files = []
            for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
                for file in files:
                    if file.endswith(".fastq.gz"):
                        fastq_gz_files.append(os.path.join(path, file))

        for file in fastq_gz_files:
            if (wildcards.tissue_name in file) and (f"_{wildcards.tag}" in file) and (f"_{wildcards.PE_SE}" in file):
                return file
    def get_trim_threads(wildcards):
        """
        Trim galore uses 9 threads on single-ended data, and 15 cores for paired-end data.
        Note: The actual trim_galore call below does not request the maximum threads given.
        Trim galore's MAN page states that it can use UP TO this many threads, however
        """
        threads = 1
        if str(wildcards.PE_SE) == "1": threads = 16
        elif str(wildcards.PE_SE) == "2": threads = 1
        elif str(wildcards.PE_SE) == "S": threads = 9
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
        input: get_trim_input
        output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            direction="{PE_SE}"
        threads: get_trim_threads
        conda: "envs/trim.yaml"
        resources:
            mem_mb = 10240,  # 10 GB
            runtime=get_trim_runtime
        shell:
            """
            # Only process on forward reads
            input_directory="$(dirname {input})"
            output_directory="$(dirname {output})"
            
            
            if [ "{params.direction}" == "1" ]; then
                file_in_1="{input}"                                                                 # Input file 1
                file_in_2="$input_directory/{params.tissue_name}_{params.tag}_2.fastq.gz"         # Input file 2
                trim_galore --paired --cores 4 -o "$(dirname {output})" "$file_in_1" "$file_in_2"
                       
                file_out_1="$output_directory/{params.tissue_name}_{params.tag}_1_val_1.fq.gz"    # final output paired end, forward read
                file_out_2="$output_directory/{params.tissue_name}_{params.tag}_2_val_2.fq.gz"    # final output paired end, reverse read
                file_rename_1="$output_directory/trimmed_{params.tissue_name}_{params.tag}_1.fastq.gz"    # final renamed output paired end, forward read
                file_rename_2="$output_directory/trimmed_{params.tissue_name}_{params.tag}_2.fastq.gz"    # final renamed output paired end, reverse read
                
                mv "$file_out_1" "$file_rename_1"
                mv "$file_out_2" "$file_rename_2"
            
            # Skip over reverse-reads. Create the output file so snakemake does not complain about the rule not generating output
            elif [ "{params.direction}" == "2" ]; then
                touch {output}
            
            # Work on single-end reads
            elif [ "{params.direction}" == "S" ]; then
                trim_galore --cores 4 -o "$(dirname {output})" "{input}"
                
                file_name="$output_directory/{params.tissue_name}_{params.tag}_S_trimmed.fq.gz"   # final output single end
                
                mv "$file_name" "{output}"
            fi
            """

    rule fastqc_trim:
        input: rules.trim.output
        output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
        params:
            file_two_input=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "trimmed_reads", "trimmed_{tissue_name}_{tag}_2.fastq.gz"),
            file_two_out=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "fastqc", "trimmed_reads", "trimmed_{tissue_name}_{tag}_2_fastqc.zip"),
            direction="{PE_SE}"
        threads: get_fastqc_threads
        conda: "envs/fastqc.yaml"
        resources:
            # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
            mem_mb=15360, # 15 GB
            runtime=get_fastqc_runtime
        shell:
            """
            output_directory="$(dirname {output})"
            mkdir -p "$output_directory"
            if [ "{params.direction}" == "1" ]; then
                fastqc {input} --threads {threads} -o "$output_directory"
                printf "FastQC finished $(basename {input}) (1/2)\n\n"
                fastqc {params.file_two_input} --threads {threads} -o "$output_directory"
                printf "FastQC finished $(basename {params.file_two_input}) (2/2)\n\n"
            elif [ "{params.direction}" == "2" ]; then
                touch {output}
            elif [ "{params.direction}" == "S" ]; then
                fastqc {input} --threads {threads} -o "$output_directory"
                printf "FastQC finished $(basename {input}) (1/1)\n\n"
            fi
            """

def get_direction_from_name(file: str):
    file_name = os.path.basename(file)
    purge_extension = file_name.split(".")[0]
    direction = purge_extension.split("_")[-1]
    return direction
def collect_star_align_input(wildcards):
    if perform_trim():
        # Have not expanded output from rule trim, need to expand it here
        in_files = sorted(expand(rules.trim.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE()))
    else:
        # already expanding output from dump_fastq, no need to expand it here
        in_files = sorted(expand(rules.dump_fastq.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE()))

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
    return min(len(input.reads) * 60 * attempt, 10079)
rule star_align:
    input:
        reads=collect_star_align_input,
        genome_dir=rules.generate_genome.output.genome_dir,
        genome_file=rules.generate_genome.output.genome_file,
        rule_complete=os.path.join(config["ROOTDIR"], "temp", "rule_complete", "generate_genome.complete")
    output: os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}.tab")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        star_output=os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "aligned_reads", "{tag}", "{tissue_name}_{tag}_ReadsPerGene.out.tab")
    threads: 40
    conda: "envs/star.yaml"
    resources:
        mem_mb=51200,# 50 GB
        runtime=get_star_align_runtime
    shell:
        """
        PREFIX="{config[ROOTDIR]}/data/{params.tissue_name}/aligned_reads/{params.tag}/{params.tissue_name}_{params.tag}_"
        printf "prefix is $PREFIX"
        STAR --runThreadN {threads} \
		--readFilesCommand "zcat" \
		--readFilesIn {input.reads} \
		--genomeDir {input.genome_dir} \
		--outFileNamePrefix $PREFIX \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--quantMode GeneCounts

		mv {params.star_output} {output}
        """


def multiqc_get_dump_fastq_data(wildcards):
    if perform_prefetch():
        output = expand(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw", "{tissue_name}_{tag}_{PE_SE}.fastq.gz"), zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        output = []
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                output.append(os.path.join(path, file))
    return_files = []
    for file in output:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files
def multiqc_get_fastqc_data(wildcards):
    if perform_trim():
        output_files = expand(rules.fastqc_trim.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())
    else:
        output_files = expand(rules.fastqc_dump_fastq.output, zip, tissue_name=get_tissue_name(), tag=get_tags(), PE_SE=get_PE_SE())

    return_files = []
    for file in output_files:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files
def multiqc_get_star_data(wildcards):
    return_files = []
    for file in expand(rules.star_align.output, zip, tissue_name=get_tissue_name(), tag=get_tags()):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files
rule multiqc:
    input:
        fastqc_data = multiqc_get_fastqc_data,
        star_data = multiqc_get_star_data,
        dump_fastq_data = multiqc_get_dump_fastq_data
    output:
        output_file = os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "multiqc", "{tissue_name}_multiqc_report.html"),
        output_directory = directory(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "multiqc"))
    params:
        input_directory=os.path.join(config["ROOTDIR"], "data", "{tissue_name}")
    threads: 1
    conda: "envs/multiqc.yaml"
    resources:
        mem_mb=1024,  # 1 GB
        runtime=lambda wildcards, attempt: 30 * (attempt * 0.75)  # 30 minutes, don't need much more time than this if it fails
    shell:
        """
        mkdir -p "{output}"
        multiqc "{params.input_directory}" --filename {wildcards.tissue_name}_multiqc_report.html --outdir {output.output_directory}
        """
