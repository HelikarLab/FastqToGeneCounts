import lib.CheckPerformRules as CheckPerformRules
import os

# def perform_prefetch():
#     if str(config["PERFORM_PREFETCH"]).lower() == "true":
#         return True
#     else:
#         return False

if CheckPerformRules.perform_prefetch():
    rule distribute_init_files:
        input: ancient(config["MASTER_CONTROL"])
        output: os.path.join(config["ROOTDIR"],"controls","init_files","{tissue_name}_{tag}.csv")
        params: id="{tissue_name}_{tag}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 1500 * attempt,
            runtime=lambda wildcards, attempt: 5 * attempt
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
        output: os.path.join(config["ROOTDIR"],"temp","prefetch","{tissue_name}_{tag}","{srr_code}.sra")
        conda: "../envs/SRAtools.yaml"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 10000 * attempt,
            runtime=lambda wildcards, attempt: 30 * attempt
        shell:
            """
            IFS=","
            while read srr name endtype; do
                # prefetch has a default max size of 20G. Effectively remove this size by allowing files up to 1TB to be downloaded
                prefetch $srr --max-size 1024000000000 --output-file {output} || touch {output}
            done < {input}
            """
