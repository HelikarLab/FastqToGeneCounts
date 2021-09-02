"""
This is an attempt to convert HCC_Data/align_liver_PE.bash into a snakefile
"""

import glob
import os
import subprocess
import csv
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
    with open(config["MASTER_INIT"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            id = line[1].split("_")[0]  # naiveB_S1R1 -> naiveB
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
    with open(config["MASTER_INIT"], "r") as rfile:
        reader = csv.reader(rfile)
        for line in reader:
            tag = line[1].split("_")[-1]
            tag_data.append(tag)  # example: S1R1

    return tag_data


rule master:
    input:
        # Distribute and download fastq
        expand(os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw"),tissue_name=get_tissue_name())



rule distribute_init_files:
    input: config["MASTER_CONTROL"]
    output: temp(os.path.join(config["ROOTDIR"], "controls", "init_files", "{tissue_name}_{tag}.csv"))
    params:
        id = "{tissue_name}_{tag}"
    # output: os.path.join(config["ROOTDIR"], "controls", "init_files", "{tissue_name}_{tag}.csv")
    run:
        # Get lines in master control file
        # Open output for writing
        lines = open(str(input), "r").readlines()
        wfile = open(str(output), "w")

        for line in lines:
            line = line.rstrip()  # remove trailing newline

            # Only write line if the output file has the current tissue-name_tag (naiveB_S1R1) in the file name
            if params.id in line:
                wfile.write(line)

        wfile.close()

rule download_fastq:
    input: expand(rules.distribute_init_files.output, tag=get_tag_data(), allow_missing=True)
    output: directory(os.path.join(config["ROOTDIR"], "data", "{tissue_name}", "raw"))
    shell:
        """
        module load SRAtoolkit/2.10 
        READ_FILE={input}
        ROOTDIR={config[ROOTDIR]}
        DATAEXT="/data/raw/"
        DATADIR="$ROOTDIR$DATAEXT"
        FQEXT=".fastq.gz
        FWD="_1"
        REV="_2"
        OLDIFS=$IFS
        IFS=","
        [ ! -f $READ_FILE ] &while read srr name endtype; do
            echo "----------------------------"
            endtype='echo "$endtype" | xargs'
            if [ $endtype == "PE" ]; the
                sfile1="$DATADIR$srr$FWD$FQEXT"
                sfile2="$DATADIR$srr$REV$FQEXT"
                rfile1="$DATADIR$name$FWD$FQEXT"
                rfile2="$DATADIR$name$REV$FQEXT"
                if [ ! -f "$rfile1" ] && [ ! -f "$rfile2" ]; then
                    fastq-dump --split-files $srr --outdir $DATADIR --gzip
                    mv $sfile1 $rfile1
                    mv $sfile2 $rfile2
                else
                    echo "already downloaded"
                fi
            elif [ $endtype == "SE" ]; then
                sfile="$DATADIR$srr$FQEXT" 
                rfile="$DATADIR$name$FQEXT"
                if [ ! -f "$rfile" ]; then
                    fastq-dump $srr --outdir $DATADIR --gzip
                    mv $sfile $rfile
                else
                    echo "already downloaded"
                fi
            else
                echo "layout not specified with 'PE' or 'SE', skipping over...."
            fi
        done < $READ_FILE
        IFS=$OLDIFS
        """
