"""
This file is responsible for holding the "perform_*" functions
this will (hopefully) help clean up the master snakefile, as it is quite lengthly with these function in place

This is required after splitting the master Snakefile into multiple sub-workflows, as these new sub-workflows
    do not have access to the functions implemented by the master Snakefile (such as def perform_prefetch)

As a result of this, we must implement these functions into a file that can be accessible by any Snakefile

"""
from lib.ConfigParser import ConfigParser
import csv
import os
from snakemake.io import glob_wildcards, expand
import warnings
import sys

config = ConfigParser("snakemake_config.yaml")


def hello_world():
    print("Hello world!")

def perform_prefetch() -> bool:
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return True
    else:
        return False


def perform_trim():  # QC
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return True
    else:
        return False


def perform_screen():  # QC
    if str(config["PERFORM_SCREEN"]).lower() == "true":
        return True
    else:
        return False


def perform_get_insert_size():
    if str(config["PERFORM_GET_INSERT_SIZE"]).lower() == "true":
        return True
    else:
        return False


def perform_get_fragment_size():  # for zFPKM QC
    if str(config["PERFORM_GET_FRAGMENT_SIZE"]).lower() == "true":
        return True
    else:
        return False


def perform_get_rnaseq_metrics():  # QC
    if str(config["PERFORM_GET_RNASEQ_METRICS"]).lower() == "true":
        return True
    else:
        return False


def get_from_master_config(attribute: str) -> list[str]:
    valid_inputs = ["SRR", "tissue", "tag", "PE_SE"]
    sub_list = ["tissue", "tag"]
    if attribute not in valid_inputs:
        sys.exit(
            f"\nInvalid attribute input. '{attribute}' is not one of: {valid_inputs}\n"
        )
    else:
        collect_attributes = []
        index_value = valid_inputs.index(attribute)

        # We have to subtract one because "tissue" and "tag" are in the same index, thus the index value in valid_inputs is increased by one
        if index_value >= 2:
            index_value -= 1
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
        fastq_input = glob_wildcards(
            os.path.join(
                config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            )
        )
        return fastq_input.tissue_name


def get_tags() -> list[str]:
    if perform_prefetch():
        return get_from_master_config("tag")
    else:
        fastq_input = glob_wildcards(
            os.path.join(
                config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            )
        )
        return fastq_input.tag


def get_PE_SE() -> list[str]:
    if perform_prefetch():
        return get_from_master_config("PE_SE")
    else:
        fastq_input = glob_wildcards(
            os.path.join(
                config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            )
        )
        return fastq_input.PE_SE


def get_sample() -> list[str]:
    if perform_prefetch():
        tag = get_from_master_config("tag")
    else:
        fastq_input = glob_wildcards(
            os.path.join(
                config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            )
        )
        tag = fastq_input.tag
    sample = []
    for t in tag:
        sample.append(t.split("R")[0])
    return sample


def perform_screen_rule(wildcards):
    """
    If screening for contamination, return fastq_screen output
    """
    if perform_screen():
        return expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "fq_screen",
                "{tissue_name}_{tag}_{PE_SE}_screen.txt",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    else:
        return []


def perform_get_insert_size_rule(wildcards):
    """
    If getting insert sizes with picard, return GetinsertSizeMetrics output
    """
    if perform_get_insert_size():
        return expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "picard",
                "insert",
                "{tissue_name}_{tag}_insert_size.txt",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
        )
    else:
        return []


def perform_get_fragment_size_rule(wildcards):
    """
    If getting fragment sizes with deeptools, return RNA_fragment_size.py output
    """
    if perform_get_fragment_size():
        return expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "fragmentSizes",
                "{tissue_name}_{tag}_fragment_length.txt",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
        )
    else:
        return []


def perform_trim_rule(wildcards):
    """
    If we are performing trimming, return trim's output
    """
    if perform_trim():
        return expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "trimmed_reads",
                "trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    else:
        return []


def fastqc_trimmed_reads(wildcards):
    """
    If we are going to trim, return output for rule fastqc_trim
    """
    if perform_trim():
        return expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "fastqc",
                "trimmed_reads",
                "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    else:
        return []


def perform_dump_fastq(wildcards):
    if perform_prefetch():
        data = expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "raw",
                "{tissue_name}_{tag}_{PE_SE}.fastq.gz",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
        return data
    else:
        return []


def get_tag(file_path: str) -> str:
    file_name = os.path.basename(file_path)
    purge_extension = file_name.split(".")[0]
    tag = purge_extension.split("_")[-1]
    return str(tag)


def get_fastqc_threads(wildcards, input):
    threads = 1
    tag = get_tag(str(input))
    if tag in ["1", "S"]:
        threads = 4
    elif tag == "2":
        threads = 2
    return threads


def get_fastqc_runtime(wildcards, input, attempt):
    tag = get_tag(str(input[0]))
    runtime = 1
    if tag in ["1", "S"]:
        runtime = 150 * attempt  # 2.5 hours
    elif tag == "2":
        runtime = 150 * attempt
    return runtime


def get_screen_input(wildcards):
    """
    aggregate filesnames of all fastqs
    """
    if perform_prefetch():
        return checkpoints.dump_fastq.get(**wildcards).output
    else:
        fastq_gz_files = []
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                if (
                        (wildcards.tissue_name in file)
                        and (f"_{wildcards.tag}" in file)
                        and (f"_{wildcards.PE_SE}" in file)
                ):
                    return os.path.join(path, file)


def get_screen_runtime(wildcards, attempt):
    """
    runtime should be relatively short since only a fraction of reads are used
    """
    # Time in minutes
    runtime = 30 * attempt
    return runtime


def get_trim_input(wildcards):
    if perform_prefetch():
        return checkpoints.dump_fastq.get(**wildcards).output
    else:
        fastq_gz_files = []
        for path, subdir, files in os.walk(config["DUMP_FASTQ_FILES"]):
            for file in files:
                if (
                        (wildcards.tissue_name in file)
                        and (f"_{wildcards.tag}" in file)
                        and (f"_{wildcards.PE_SE}" in file)
                ):
                    return os.path.join(path, file)


def get_trim_threads(wildcards):
    """
    Trim galore uses 9 threads on single-ended data, and 15 cores for paired-end data.
    Note: The actual trim_galore call below does not request the maximum threads given.
    Trim galore's MAN page states that it can use UP TO this many threads, however
    """
    threads = 1
    if str(wildcards.PE_SE) == "1":
        threads = 16
    elif str(wildcards.PE_SE) == "2":
        threads = 1
    elif str(wildcards.PE_SE) == "S":
        threads = 9
    return threads


def get_trim_runtime(wildcards, attempt):
    """
    Trim galore takes more time on single-ended data and the forward read of paired end data
    It only touches the output file of the reverse-read paired-end data
    """
    runtime = 5  # minutes
    if str(wildcards.PE_SE) == "1":
        runtime = 120 * attempt
    elif str(wildcards.PE_SE) == "2":
        runtime = 120 * attempt
    elif str(wildcards.PE_SE) == "S":
        runtime = 120 * attempt
    return runtime


def get_direction_from_name(file: str):
    file_name = os.path.basename(file)
    purge_extension = file_name.split(".")[0]
    direction = purge_extension.split("_")[-1]
    return direction


def collect_star_align_input(wildcards):
    if perform_trim():
        # Have not expanded output from rule trim, need to expand it here
        in_files = sorted(
            expand(
                rules.trim.output,
                zip,
                tissue_name=get_tissue_name(),
                tag=get_tags(),
                PE_SE=get_PE_SE(),
            )
        )
    else:
        # already expanding output from dump_fastq, no need to expand it here
        # in_files = sorted(expand(rules.dump_fastq.output,zip,tissue_name=get_tissue_name(),tag=get_tags(),PE_SE=get_PE_SE()))
        in_files = sorted(
            expand(
                checkpoints.dump_fastq.output,
                zip,
                tissue_name=get_tissue_name(),
                tag=get_tags(),
                PE_SE=get_PE_SE(),
            )
        )

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
            elif direction == "2":
                continue
            else:
                warnings.warn(
                    f"{in_file} expects additional paired-end read! Skipping...."
                )
                continue

        if direction == "S":
            grouped_reads.append(in_file)
        elif direction == "1" and next_direction == "2":
            if (
                    in_file[:-10] == next_file[:-10]
            ):  # remove _1.fastq.gz to make sure they are same replicate
                both_reads = " ".join([in_file, next_file])
                grouped_reads.append(both_reads)
            else:
                warnings.warn(
                    f"{in_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {in_file} and one end of {next_file} failed to download"
                )

        elif direction == "1" and not next_direction == "2":
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
    This is much like what was done in the function get_dump_fastq_output
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
    return min(len(input.reads) * 60 * 4 * attempt, 10079)


def insert_size_get_star_data(wildcards):
    return_files = []
    for file in expand(
            rules.star_align.output.bam_file,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_dump_fastq_data(wildcards):
    if perform_prefetch():
        output = expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "raw",
                "{tissue_name}_{tag}_{PE_SE}.fastq.gz",
            ),
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
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
        output_files = expand(
            rules.fastqc_trim.output,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    else:
        output_files = expand(
            rules.fastqc_dump_fastq.output,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    return_files = []
    for file in output_files:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_star_data(wildcards):
    return_files = []
    for file in expand(
            rules.star_align.output.gene_table,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_screen_data(wildcards):
    if perform_screen():
        output_files = expand(
            rules.contaminant_screen.output,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
            PE_SE=get_PE_SE(),
        )
    else:
        output_files = []
    return_files = []
    for file in output_files:
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_insertsize_data(wildcards):
    return_files = []
    if perform_get_insert_size():
        for file in expand(
                rules.get_insert_size.output.txt,
                zip,
                tissue_name=get_tissue_name(),
                tag=get_tags(),
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_fragmentsize_data(wildcards):
    return_files = []
    if perform_get_fragment_size():
        for file in expand(
                rules.get_fragment_size.output,
                zip,
                tissue_name=get_tissue_name(),
                tag=get_tags(),
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_rnaseq_data(wildcards):
    return_files = []
    for file in expand(
            rules.get_rnaseq_metrics.output,
            zip,
            tissue_name=get_tissue_name(),
            tag=get_tags(),
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files
