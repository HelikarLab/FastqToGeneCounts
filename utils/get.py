import os
import csv
import sys
import snakemake
from pathlib import Path
from . import perform
from .constants import EndType
from typing import Union

def from_master_config(config: dict, attribute: str) -> list[str]:
    valid_inputs = ["SRR", "tissue", "tag", "PE_SE"]
    sub_list = ["tissue", "tag"]
    if attribute not in valid_inputs:
        sys.exit(f"\nInvalid attribute input. '{attribute}' is not one of: {valid_inputs}\n")
    else:
        collect_attributes = []
        index_value = valid_inputs.index(attribute)

        # We have to subtract one because "tissue" and "tag" are in the same index, thus the index value in valid_inputs is increased by one
        if index_value >= 2:
            index_value -= 1

        control_lines = open(config["MASTER_CONTROL"], "r").readlines()
        reader = csv.reader(control_lines)

        for line in reader:

            # Get the column from master_control we are interested in
            column_value = line[index_value]
            PE_SE_value = EndType[line[2]]  # PE, SE, or SLC

            # test if we are looking for "tissue" or "tag", as these two values are located at master_control index 1
            if attribute in sub_list:
                sub_index = sub_list.index(attribute)
                split_list = str(line[index_value]).split("_")

                # We must append the target attribute twice if it is paired end, once if it is single end
                if PE_SE_value in [EndType.PE, EndType.SLC]:
                    target_attribute = [split_list[sub_index], split_list[sub_index]]
                elif PE_SE_value == EndType.SE:
                    target_attribute = [split_list[sub_index]]

            # Test if we are gathering the ended-ness (PE, SE, SLC)
            elif attribute == "PE_SE":
                # We must append the target attribute twice if it is paired end, once if it is single end
                if column_value in [EndType.PE.name, EndType.SLC.name]:  # paired end or single cell
                    target_attribute = ["1", "2"]
                elif column_value == EndType.SE.name:  # Single end
                    target_attribute = ["S"]

            # If we are doing anything else, simply append the column value the appropriate number of times
            else:
                if PE_SE_value in [EndType.PE, EndType.SLC]:
                    target_attribute = [line[index_value], line[index_value]]
                elif PE_SE_value == EndType.SE:
                    target_attribute = [line[index_value]]

            collect_attributes += target_attribute

        return collect_attributes


def srr_code(config: dict) -> list[str]:
    """
    Only should be getting SRR values if we are performing prefetch
    """
    if perform.trim(config=config):
        return from_master_config(config=config, attribute="SRR")


def tissue_name(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="tissue")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tissue_name


def tags(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="tag")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.tag


def PE_SE(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        return from_master_config(config=config, attribute="PE_SE")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        return fastq_input.PE_SE

def tag_from_filename(file_path: Union[str, Path]) -> str:
    file_name = os.path.basename(file_path)
    purge_extension = file_name.split(".")[0]
    tag = purge_extension.split("_")[-1]
    return str(tag)

def direction_from_name(file: str):
    file_name = os.path.basename(file)
    purge_extension = file_name.split(".")[0]
    direction = purge_extension.split("_")[-1]
    return direction

def sample(config: dict) -> list[str]:
    if perform.prefetch(config=config):
        tag = from_master_config(config=config, attribute="tag")
    else:
        fastq_input = snakemake.io.glob_wildcards(os.path.join(config["DUMP_FASTQ_FILES"], "{tissue_name}_{tag}_{PE_SE}.fastq.gz"))
        tag = fastq_input.tag

    sample = []
    for t in tag:
        sample.append(t.split("R")[0])
    return sample
