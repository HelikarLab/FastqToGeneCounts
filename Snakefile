"""
This branch is an attempt to move portions of the workflow into "sub-workflows"
for easier understanding of what is happening throuhgout this pipeline
"""
from lib.PerformFunctions import *
import sys

configfile: "snakemake_config.yaml"

# ---- Include sub-workflows we need ----
include: "workflows/preliminary_setup.smk"
include: "workflows/MADRID_copy_setup.smk"

# Validate users are using conda. This is important for temporary conda environments defined in the workflow
# From: https://stackoverflow.com/a/71678773/13885200
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass the '--use-conda' flag to snakemake.\nExample: snakemake --cores 10 --use-conda\n\n")
    sys.exit(1)


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


rule_all = [
    "preroundup.txt",  # pre-roundup
    config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],  # Generate Genome
    perform_dump_fastq,  # dump_fastq
    perform_screen_rule,  # fastq_screen
    perform_trim_rule,  # trim reads
    fastqc_trimmed_reads,

    # FastQC
    expand(
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags(),
        PE_SE=get_PE_SE()
    ),


    # # STAR aligner
    expand(
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "aligned_reads",
            "{tag}",
            "{tissue_name}_{tag}.tab"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags()
    ),

    expand(
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "aligned_reads",
            "{tag}",
            "{tissue_name}_{tag}.bam.bai"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags()
    ),

    # copy .tab
    expand(
        os.path.join(
            "MADRID_input",
            "{tissue_name}",
            "geneCounts",
            "{sample}",
            "{tissue_name}_{tag}.tab"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags(),
        sample=get_sample()
    ),

    # get rnaseq metrics
    expand(
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "picard",
            "rnaseq",
            "{tissue_name}_{tag}_rnaseq.txt"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags()
    ),

    # copy strandedness
    expand(
        os.path.join(
            "MADRID_input",
            "{tissue_name}",
            "strandedness",
            "{sample}",
            "{tissue_name}_{tag}_strandedness.txt"),
        zip,
        tissue_name=get_tissue_name(),
        tag=get_tags(),
        sample=get_sample()
    ),

    # MultiQC
    expand(
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "multiqc",
            "{tissue_name}_multiqc_report.html"),
        tissue_name=get_tissue_name()
    )
]

if perform_get_insert_size():
    rule_all.extend([
        # Get Insert sizes
        perform_get_insert_size_rule,

        # copy insert
        expand(
            os.path.join(
                "MADRID_input",
                "{tissue_name}",
                "insertSizeMetrics",
                "{sample}",
                "{tissue_name}_{tag}_insert_size.txt"),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            sample=get_sample()
        )
    ])

if perform_get_fragment_size():
    rule_all.extend([
        # get fragment lengths
        perform_get_fragment_size_rule,

        # copy fragment
        expand(
            os.path.join(
                "MADRID_input",
                "{tissue_name}",
                "fragmentSizes",
                "{sample}",
                "{tissue_name}_{tag}_fragment_size.txt"),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            sample=get_sample()
        )
    ])

rule all:
    input: rule_all


def fastqc_dump_fastq_input(wildcards):
    """
    This function will return the input for fastqc_dump_fastq
    It is going to return forward read AND reverse read if the input file is the forward read
    If input is the reverse read, it will only return the reverse read
    If input is a single end read, it will only return the single end read
    """
    if perform_prefetch():
        checkpoint_output = str(checkpoints.dump_fastq.get(**wildcards).output)
        if str(wildcards.PE_SE) == "1":
            file_two = os.path.join(os.path.dirname(checkpoint_output),f"{wildcards.tissue_name}_{wildcards.tag}_2.fastq.gz")
            return [checkpoint_output, file_two]
        else:
            return checkpoint_output
    else:
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                if (wildcards.tissue_name in file) and (wildcards.tag in file) and (f"_{wildcards.PE_SE}" in file):
                    file_one = os.path.join(path,file)
                    if str(wildcards.PE_SE) == "1":
                        file_two = os.path.join(path,f"{wildcards.tissue_name}_{wildcards.tag}_2.fastq.gz")
                        return [file_one, file_two]
                    else:
                        return file_one

rule fastqc_dump_fastq:
    input: fastqc_dump_fastq_input
    output:
        os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
        ),
    params:
        file_one_zip=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
        ),
        file_one_html=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "{tissue_name}_{tag}_{PE_SE}_fastqc.html",
        ),
        file_two_zip=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "{tissue_name}_{tag}_2_fastqc.zip",
        ),
        file_two_html=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "{tissue_name}_{tag}_2_fastqc.html",
        ),
        file_one_zip_rename=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
        ),
        file_one_html_rename=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.html",
        ),
        file_two_zip_rename=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_2_fastqc.zip",
        ),
        file_two_html_rename=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_2_fastqc.html",
        ),
    threads: get_fastqc_threads
    conda: "envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 15000 * attempt,# 15 GB
        runtime=get_fastqc_runtime  # 150 minutes * attempt number
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory" || \
            ( touch {params.file_one_zip} && touch {params.file_one_html} \
              touch {params.file_two_zip} && touch {params.file_two_html} )

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
            mv "{params.file_two_zip}" "{params.file_two_zip_rename}"
            mv "{params.file_two_html}" "{params.file_two_html_rename}"

        elif [ "{wildcards.PE_SE}" == "2" ]; then
            touch {output}

        elif [ "{wildcards.PE_SE}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory" || \
            ( touch {params.file_one_zip} && touch {params.file_one_html} )

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
        fi
        """


if perform_screen():
    rule contaminant_screen:
        input:
            files=get_screen_input,
            genomes=rules.get_screen_genomes.output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fq_screen","{tissue_name}_{tag}_{PE_SE}_screen.txt")
        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            direction="{PE_SE}"
        conda: "envs/screen.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 5000 * attempt,# 5 GB
            runtime=get_screen_runtime
        shell:
            """
            fastq_screen --aligner Bowtie2 --conf FastQ_Screen_Genomes/fastq_screen.conf {input} || touch {params.tissue_name}_{params.tag}_{params.direction}_screen.txt
            mv {params.tissue_name}_{params.tag}_{params.direction}_screen.txt results/data/{params.tissue_name}/fq_screen/
            """

if perform_trim():
    rule trim:
        input: get_trim_input
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
        params:
            tissue_name="{tissue_name}",
            tag="{tag}",
            direction="{PE_SE}"
        threads: get_trim_threads
        conda: "envs/trim.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 10000 * attempt,# 10 GB
            runtime=get_trim_runtime
        shell:
            """
            # Only process on forward reads
            input_directory="$(dirname {input})"
            output_directory="$(dirname {output})"

            if [ "{params.direction}" == "1" ]; then
                file_in_1="{input}"                                                               # Input file 1
                file_in_2="$input_directory/{params.tissue_name}_{params.tag}_2.fastq.gz"         # Input file 2

                file_out_1="$output_directory/{params.tissue_name}_{params.tag}_1_val_1.fq.gz"    # final output paired end, forward read
                file_out_2="$output_directory/{params.tissue_name}_{params.tag}_2_val_2.fq.gz"    # final output paired end, reverse read
                file_rename_1="$output_directory/trimmed_{params.tissue_name}_{params.tag}_1.fastq.gz"    # final renamed output paired end, forward read
                file_rename_2="$output_directory/trimmed_{params.tissue_name}_{params.tag}_2.fastq.gz"    # final renamed output paired end, reverse read

                trim_galore --paired --cores 4 -o "$(dirname {output})" "$file_in_1" "$file_in_2" || (touch $file_out_1 & touch $file_out_2)

                mv "$file_out_1" "$file_rename_1"
                mv "$file_out_2" "$file_rename_2"

            # Skip over reverse-reads. Create the output file so snakemake does not complain about the rule not generating output
            elif [ "{params.direction}" == "2" ]; then
                touch {output}

            # Work on single-end reads
            elif [ "{params.direction}" == "S" ]; then
                file_out="$output_directory/{params.tissue_name}_{params.tag}_S_trimmed.fq.gz"   # final output single end
                #file_rename="$output_directory/trimmed_{params.tissue_name}_{params.tag}_S.fastq.gz"

                trim_galore --cores 4 -o "$(dirname {output})" "{input}" || touch $file_out

                mv "$file_out" "{output}"
            fi
            """


    rule fastqc_trim:
        input: rules.trim.output
        output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
        params:
            file_two_input=os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_2.fastq.gz"),
            file_two_out=os.path.join(config[
                "ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_2_fastqc.zip")
        threads: get_fastqc_threads
        conda: "envs/fastqc.yaml"
        resources:
            # fastqc allocates 250MB per thread. 250*5 = 1250MB ~= 2GB for overhead
            mem_mb=lambda wildcards, attempt, threads: attempt * threads * 2,# 15 GB
            runtime=get_fastqc_runtime
        shell:
            """
            output_directory="$(dirname {output})"
            mkdir -p "$output_directory"

            if [ "{wildcards.PE_SE}" == "1" ]; then
                fastqc {input} --threads {threads} -o "$output_directory" || touch {output}
                printf "FastQC finished $(basename {input}) (1/2)\n\n"

                fastqc {params.file_two_input} --threads {threads} -o "$output_directory" || touch {params.file_two_out}
                printf "FastQC finished $(basename {params.file_two_input}) (2/2)\n\n"

            elif [ "{wildcards.PE_SE}" == "2" ]; then
                touch {output}

            elif [ "{wildcards.PE_SE}" == "S" ]; then
                fastqc {input} --threads {threads} -o "$output_directory" || touch {output}
                printf "FastQC finished $(basename {input}) (1/1)\n\n"
            fi
            """


rule star_align:
    input:
        reads=collect_star_align_input,
        genome_dir=rules.generate_genome.output.genome_dir,
        generate_genome_complete=rules.generate_genome.output.rule_complete
    output:
        gene_table=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab"),
        bam_file=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        gene_table_output=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_ReadsPerGene.out.tab"),
        bam_output=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_Aligned.sortedByCoord.out.bam")
    threads: 40
    conda: "envs/star.yaml"
    resources:
        mem_mb=50000,# 50 GB
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
		--quantMode GeneCounts || touch {params.gene_table_output}
		mv {params.gene_table_output} {output.gene_table}
		mv {params.bam_output} {output.bam_file}
        """

rule index_bam_file:
    input: rules.star_align.output.bam_file
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt,# 1 GB
        runtime=1
    conda: "envs/samtools.yaml"
    shell:
        """
        samtools index {input} {output} || touch {output}
        """


rule get_rnaseq_metrics:
    input:
        bam=rules.star_align.output.bam_file,
        tab=rules.star_align.output.gene_table
    output:
        metrics=os.path.join(config["ROOTDIR"],"data","{tissue_name}","picard","rnaseq","{tissue_name}_{tag}_rnaseq.txt"),
        strand=os.path.join(config["ROOTDIR"],"data","{tissue_name}","strand","{tissue_name}_{tag}_strand.txt")
    params:
        ref_flat=config["REF_FLAT_FILE"],
        ribo_int_list=config["RRNA_INTERVAL_LIST"]
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * 5 * attempt,# 5 GB / attempt
        runtime=lambda wildcards, attempt: 60 * attempt
    conda: "envs/picard.yaml"
    shell:
        """
        colsums=$(grep -v "N_" {input.tab} | awk '{{unst+=$2;forw+=$3;rev+=$4}}END{{print unst,forw,rev}}') || colsums="0 1 2"
        IFS=" "
        read -ra arr <<< "$colsums"

        cnt=0
        # all values are at least one to prevent divide by zero
        unst=1
        fwd=1
        rev=1
        declare -i unst
        declare -i fwd
        declare -i rev

        for val in "${{arr[@]}}"; do
        cnt=$((cnt+=1))
        declare -i val
        if [[ $cnt == 1 ]]; then
            unst+=$val
        elif [[ $cnt == 2 ]]; then
            fwd+=$val
        elif [[ $cnt == 3 ]]; then
            rev+=$val
        fi
        done

        if [[ $(( rev / (fwd+1) )) -gt 2 ]]; then
            str_spec="SECOND_READ_TRANSCRIPTION_STRAND"
        elif [[ $(( fwd / (rev+1) )) -gt 2 ]]; then
            str_spec="FIRST_READ_TRANSCRIPTION_STRAND"
        else
            str_spec="NONE"
        fi
        echo $str_spec > {output.strand}

        picard CollectRnaSeqMetrics I={input.bam} O={output.metrics} REF_FLAT={params.ref_flat} STRAND_SPECIFICITY=$str_spec RIBOSOMAL_INTERVALS={params.ribo_int_list} || touch {output}
        """

if perform_get_insert_size():
    rule get_insert_size:
        input:
            bam=rules.star_align.output.bam_file,
            preround=rules.preroundup.output
        output:
            txt=os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","picard","insert","{tissue_name}_{tag}_insert_size.txt"),
            pdf=os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","picard","hist","{tissue_name}_{tag}_insert_size_histo.pdf")
        params:
            layout=os.path.join(config["ROOTDIR"],"data","{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt")
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: 1000 * 5 * attempt,# 5 GB / attempt
            runtime=lambda wildcards, attempt: 60 * attempt
        conda: "envs/picard.yaml"
        #wrapper:
        #    "v1.0.0/bio/picard/collectinsertsizemetrics"
        shell:
            """
            lay=$(cat {params.layout})
            if [ $lay == "paired-end"]; then
                picard CollectinsertSizeMetrics \
                I={input.bam} \
                O={output.txt} \
                H={output.pdf} \
                M=0.05 || picard CollectinsertSizeMetrics I={input} O={output.txt} H={output.pdf} M=0.5
            else
                echo "cannot collect metrics for single-end data" > {output.txt}
                touch {output.pdf}
            fi
            """

if perform_get_fragment_size():
    rule get_fragment_size:
        input:
            bam=rules.star_align.output.bam_file,
            bai=rules.index_bam_file.output
        output:
            os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","fragmentSizes","{tissue_name}_{tag}_fragment_length.txt")
        params:
            layout=os.path.join(config["ROOTDIR"],"data","{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
            bed=config["BED_FILE"]
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: 1000 * 5 * attempt,# 5 GB / attempt
            runtime=lambda wildcards, attempt: 60 * attempt
        conda: "envs/rseqc.yaml"
        shell:
            """
            files=(.snakemake/conda/*/bin/RNA_fragment_size.py) # get matches of script file ( should only be one but )      
            python3 ${{files[0]}} -r {params.bed} -i {input.bam} > {output} || touch {output}  # run first match
            """


rule multiqc:
    input:
        fastqc_data=multiqc_get_fastqc_data,
        star_data=multiqc_get_star_data,
        dump_fastq_data=multiqc_get_dump_fastq_data,
        screen_data=multiqc_get_screen_data,
        insertsize_data=multiqc_get_insertsize_data,
        rnaseq_data=multiqc_get_rnaseq_data,
        fragment_size_data=multiqc_get_fragmentsize_data
    output:
        output_file=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","multiqc","{tissue_name}_multiqc_report.html"),
        output_directory=directory(os.path.join(config["ROOTDIR"],"data","{tissue_name}","multiqc"))
    params:
        input_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}")
    threads: 1
    conda: "envs/multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt,  # 1 GB / attempt
        runtime=lambda wildcards, attempt: int(30 * (attempt * 0.75))  # 30 minutes, don't need much more time than this if it fails
    shell:
        """
        mkdir -p "{output}"
        multiqc "{params.input_directory}" --filename {wildcards.tissue_name}_multiqc_report.html --outdir {output.output_directory}
        if ls ./*.txt 1> /dev/null 2>&1; then
            rm *.txt
        fi
        if ls ./*.html 1> /dev/null 2>&1; then
            rm *.html
        fi
        if ls ./*.png 1> /dev/null 2>&1; then
            rm *.png
        fi
        if ls ./*.fastq 1> /dev/null 2>&1; then
            rm *.fastq
        fi
        """

