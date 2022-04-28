# from lib.PerformFunctions import *
import os

include: "../lib/PerformFunctions.smk"

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
    """
    # Max time is 10,080 minutes (7 days), do not let this function return more than that amount of time
    return 45 * attempt


def dump_fastq_input(wildcards):
    output_files = expand(
        rules.prefetch.output,
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags(),
        srr_code=get_srr_code(),
    )
    for file in output_files:
        if (wildcards.tissue_name in file) and (wildcards.tag in file):
            return file


def get_dump_fastq_threads(wildcards):
    """Get threads for dump fastq"""
    threads = 1
    if str(wildcards.PE_SE) in ["1", "S"]:
        threads = 40
    elif str(wildcards.PE_SE) == "2":
        threads = 1
    return threads


def get_dump_fastq_srr_code(wildcards, input):
    """Get SRR codes corresponding to dump_fastq output"""
    file_name = os.path.basename(str(input))
    srr_code = file_name.split(".")[0]
    return srr_code


if perform_prefetch():

    rule distribute_init_files:
        input:
            ancient(config["MASTER_CONTROL"]),
        output:
            os.path.join(
                config["ROOTDIR"], "controls", "init_files", "{tissue_name}_{tag}.csv"
            ),
        params:
            id="{tissue_name}_{tag}",
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 1500 * attempt,
            runtime=lambda wildcards, attempt: 5 * attempt,
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
        input:
            rules.distribute_init_files.output,
        output:
            os.path.join(
                config["ROOTDIR"],
                "temp",
                "prefetch",
                "{tissue_name}_{tag}",
                "{srr_code}.sra",
            ),
        conda:
            "../envs/SRAtools.yaml"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 10000 * attempt,
            runtime=lambda wildcards, attempt: 30 * attempt,
        shell:
            """
            IFS=","
            while read srr name endtype; do
                # prefetch has a default max size of 20G. Effectively remove this size by allowing files up to 1TB to be downloaded
                prefetch $srr --max-size 1024000000000 --output-file {output} || touch {output}
            done < {input}
            """

    checkpoint dump_fastq:
        input:
            dump_fastq_input,
        output:
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "raw",
                "{tissue_name}_{tag}_{PE_SE}.fastq.gz",
            ),
        params:
            srr_code=lambda wildcards, input: get_dump_fastq_srr_code(wildcards, input),
        threads: get_dump_fastq_threads
        conda:
            "../envs/SRAtools.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 20000 * attempt,  # 20 GB
            runtime=get_dump_fastq_runtime,
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
# TODO: convert this input to snakemake's HTTP download


# https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#read-only-web-http-s
# Not 100% sure if this is needed
rule preroundup:
    input:
        config["MASTER_CONTROL"],
    output:
        "preroundup.txt",
    params:
        rootdir=config["ROOTDIR"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 50 * attempt,
        runtime=lambda wildcards, attempt: 1 * attempt,
    shell:
        """
        IFS=","
        while read srr name endtype prep; do
            tissue=$(echo $name | cut -d '_' -f1)
            mkdir -p MADRID_input/${{tissue}}/geneCounts/
            mkdir -p MADRID_input/${{tissue}}/insertSizeMetrics/
            mkdir -p {params.rootdir}/data/${{tissue}}/layouts/
            mkdir -p {params.rootdir}/data/${{tissue}}/prepMethods/
            mkdir -p MADRID_input/${{tissue}}/layouts/
            mkdir -p MADRID_input/${{tissue}}/fragmentSizes/
            mkdir -p MADRID_input/${{tissue}}/prepMethods/
            prepl=$( echo "$prep" | tr '[:upper:]' '[:lower:]' )
            study=$( echo $name | grep -oP "_\KS\d+(?=R\d+[r]?[\d+]?)" )
            mkdir -p MADRID_input/${{tissue}}/layouts/${{study}}/
            mkdir -p MADRID_input/${{tissue}}/prepMethods/${{study}}/
            if [[ $endtype == "SE" ]]; then
                echo "single-end" > {params.rootdir}/data/${{tissue}}/layouts/${{name}}_layout.txt
                echo "single-end" > MADRID_input/${{tissue}}/layouts/${{study}}/${{name}}_layout.txt
            elif [[ $endtype == "PE" ]]; then
                echo "paired-end" > {params.rootdir}/data/${{tissue}}/layouts/${{name}}_layout.txt
                echo "paired-end" > MADRID_input/${{tissue}}/layouts/${{study}}/${{name}}_layout.txt
            else
                echo "invalid layout"
            fi
            if [[ $prepl == "mrna" ]]; then
                echo "mrna" > {params.rootdir}/data/${{tissue}}/prepMethods/${{name}}_prep_method.txt
                echo "mrna" > MADRID_input/${{tissue}}/prepMethods/${{study}}/${{name}}_prep_method.txt
            elif [[ $prep == "total" ]]; then
                echo "total" > {params.rootdir}/data/${{tissue}}/prepMethods/${{name}}_prep_method.txt
                echo "total" > MADRID_input/${{tissue}}/prepMethods/${{study}}/${{name}}_prep_method.txt
            else
                echo "invalid library preparation method. Must be total or mrna"
            fi 
        done < {input}
        touch "preroundup.txt"
        """


rule generate_genome:
    input:
        genome_fasta_file=config["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file=config["GENERATE_GENOME"]["GTF_FILE"],
    output:
        genome_dir=directory(os.path.join(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"])),
        rule_complete=touch(
            os.path.join(
                config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],
                "generate_genome.complete",
            )
        ),
    threads: 40
    resources:
        mem_mb=50000,  # 50 GB
        runtime=lambda wildcards, attempt: 120 * attempt,
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99


        """


if perform_screen():
    rule get_screen_genomes:
        """
        Download genomes to screen against
        """
        output:
            directory("FastQ_Screen_Genomes"),
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 500 * attempt,
            runtime=lambda wildcards, attempt: 240 * attempt,
        conda:
            "../envs/screen.yaml"
        shell:
            """
            if [[ ! -d "./FastQ_Screen_Genomes" ]]; then
                fastq_screen --get_genomes
                sed -i 's/\/data1\///' FastQ_Screen_Genomes/fastq_screen.conf # remove data1/ from screen genome paths
            else
                touch -c ./FastQ_Screen_Genomes/*
            fi
            """
