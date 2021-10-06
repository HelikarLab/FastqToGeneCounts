def generate_output_tuples(output_list: list[str]):
    """
    This function will generate a list of tuples that group like-files together
    Example:
        input:
            [
             "results/data/naiveB/raw/naiveB_S1R1_1.fastq.gz",
             "results/data/naiveB/raw/naiveB_S1R1_2.fastq.gz",
             "results/data/naiveB/raw/naiveB_S1R2.fastq.gz"
            ]
        output:
            [
                ("results/data/naiveB/raw/naiveB_S1R1_1.fastq.gz", "results/data/naiveB/raw/naiveB_S1R1_2.fastq.gz"),
                ("results/data/naiveB/raw/naiveB_S1R2.fastq.gz")
            ]
    :param output_list: A list of strings containing output file locations
    :return: A list of tuples containing grouped output file locations
    """

    new_list = []
    for i, output_file in enumerate(output_list):
        id = output_file.split("/")[-1].strip(".fastq.gz")
        try:  # handle final index
            next_file = output_list[i + 1]
            next_id = next_file.split("/")[-1].strip(".fastq.gz")
        except:
            if id.endswith("_S"):
                new_list.append(output_file)
                continue
            elif id.endswith("_2"):
                continue
            else:
                warnings.warn(f"{output_file} expects additional paired-end read! Skipping....")
                continue

        if id.endswith("_2"):
            continue  # skip reverse reads if not accompanied by their forwar
        elif id.endswith("_1") and next_id.endswith("_2"):
            if id[:-3] == next_id[:-3]:
                new_list.append((output_file, output_list[i + 1]))
            else:
                warnings.warn(f"{output_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {id} and one end of {next_id} failed to download")
        elif id.endswith("_S"):
            new_list.append(output_file)
        elif id.endswith("_1") and not next_id.endswith("_2"):
            warnings.warn(f"{output_file} expects additional paired-end read, it may have failed to download! Skipping....")
        else:
            warnings.warn(f"{output_file} not handled, unknown reason!")

    return new_list


import numpy as np
import subprocess
import os
import warnings

input_list = sorted(str(snakemake.input).split(" "))
snakemake_output = sorted(str(snakemake.output).split(" "))
output_list = generate_output_tuples(snakemake_output)

# Get unique items from list in the original order they were added
# input_index = np.unique(input_list, return_index=True)[1]
# input_list = [input_list[i] for i in sorted(input_index)]
# output_list = generate_output_tuples(output_list)

# iterate through input and output items
for i, (in_file, out_files) in enumerate(zip(input_list, output_list)):
    if type(out_files) is tuple:
        out_directory = os.path.dirname(out_files[0])
        os.makedirs(out_directory, exist_ok=True)
        subprocess.run(["parallel-fastq-dump", "--sra-id", str(in_file), "--threads", str(snakemake.threads), "--outdir", str(out_directory), "--gzip", "--split-files"])

        # fastq_dumped_files pulls ALL files in output directory
        # This is not good, as it also includes correctly-named output files
        # We want to include only files that have "SRR" in the name
        fastq_dumped_files = [file for file in sorted(os.listdir(out_directory)) if "SRR" in file]
        for j, (old_file, new_file) in enumerate(zip(fastq_dumped_files, out_files)):
            old_file_path = os.path.join(out_directory, old_file)
            os.rename(old_file_path, new_file)
    else:
        out_directory = os.path.dirname(out_files)
        os.makedirs(out_directory, exist_ok=True)
        subprocess.run(["parallel-fastq-dump", "--sra-id", str(in_file), "--threads", str(snakemake.threads), "--outdir", str(out_directory), "--gzip"])

        # fastq_dumped_files pulls ALL files in output directory
        # This is not good, as it also includes correctly-named output files
        # We want to include only files that have "SRR" in the name
        fastq_dumped_files = [file for file in sorted(os.listdir(out_directory)) if "SRR" in file]
        old_file_path = os.path.join(out_directory, fastq_dumped_files[0])
        os.rename(old_file_path, str(out_files))
