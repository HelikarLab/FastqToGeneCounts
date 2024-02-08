import os
import csv
import warnings
import pandas as pd
from pathlib import Path
import snakemake
from snakemake import io
from utils import get, perform, validate
from utils.constants import Layout, PrepMethod

configfile: "config.yaml"

# Validate file before reading with pandas
if validate.validate(config=config):
    print(f"Control file ({config['MASTER_CONTROL']}) is valid! Continuing...")

os.makedirs(config["ROOTDIR"],exist_ok=True)

# Get the delimiter from the master control file; from: https://stackoverflow.com/questions/16312104
dialect = csv.Sniffer().sniff(open(config["MASTER_CONTROL"]).read(1024))
samples: pd.DataFrame = pd.read_csv(
    config["MASTER_CONTROL"],
    delimiter=str(dialect.delimiter),
    names=["srr", "sample", "endtype", "prep_method"]
)

config_file_basename = os.path.basename(config["MASTER_CONTROL"]).split(".")[0]
contaminant_screen_directories = [
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Adapters"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Arabidopsis"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Drosophila"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"E_coli"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Human"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Lambda"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Mitochondria"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Mouse"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"PhiX"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Rat"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Vectors"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Worm"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"Yeast"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"rRNA"),
    os.path.join(config["CONTAMINANT_GENOME_PATH"],"fastq_screen.conf")
]


def rule_all_prefetch(wildcards):
    if not perform.prefetch(config=config):
        return []
    return expand(
        os.path.join(
            config["ROOTDIR"],"temp",
            "prefetch","{tissue_name}","{tissue_name}_{tag}",
            "{tissue_name}_{tag}.sra"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    )


def rule_all_fastq_trimmed_reads(wildcards):
    """
    If we are going to trim, return output for rule fastqc_trim
    """
    if not perform.trim(config=config):
        return []
    return expand(
        os.path.join(config["ROOTDIR"],"data",
            "{tissue_name}","fastqc","trimmed_reads",
            "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"
        ),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        PE_SE=get.PE_SE(config=config)
    )


def rule_all_contaminant_screen(wildcards):
    """
    If screening for contamination, return fastq_screen output
    """
    if not perform.screen(config=config):
        return []
    return snakemake.io.expand(
        os.path.join(config["ROOTDIR"],"data","{tissue_name}","fq_screen","{tissue_name}_{tag}_{PE_SE}_screen.txt"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        PE_SE=get.PE_SE(config=config)
    )


def rule_all_insert_sizes(wildcards):
    """
    If getting insert sizes with picard, return GetinsertSizeMetrics output
    """
    if not perform.get_insert_size(config=config):
        return []
    return snakemake.io.expand(
        os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","picard","insert","{tissue_name}_{tag}_insert_size.txt"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    )


def rule_all_fragment_sizes(wildcards):
    """
    If getting fragment sizes with deeptools, return RNA_fragment_size.py output
    """
    if not perform.get_fragment_size(config=config):
        return []
    return snakemake.io.expand(
        os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","fragmentSizes","{tissue_name}_{tag}_fragment_length.txt"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    )


def rule_all_trim(wildcards):
    """
    If we are performing trimming, return trim's output
    """
    if not perform.trim(config=config):
        return []
    return snakemake.io.expand(
        os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        PE_SE=get.PE_SE(config=config)
    )


def rule_all_dump_fastq(wildcards):
    if not perform.prefetch(config=config) or not perform.dump_fastq(config=config):
        return []
    
    return expand(
        os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config),
        PE_SE=get.PE_SE(config=config)
    )


def rule_all_screen_genomes(wildcards):
    if not perform.screen(config=config):
        return []
    return contaminant_screen_directories


rule all:
    input:
        config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],
        rule_all_screen_genomes,
        rule_all_prefetch,
        rule_all_dump_fastq,
        rule_all_trim,
        rule_all_fastq_trimmed_reads,
        expand(os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ),
        
        expand(os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","prepMethods","{tissue_name}_{tag}_prep_method.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ),
        
        expand(
            os.path.join(config[
                "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        ),
        
        expand(
            os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ),
        
        expand(
            os.path.join("COMO_input","{tissue_name}","geneCounts","{sample}","{tissue_name}_{tag}.tab"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            sample=get.sample(config=config)
        ),
        
        expand(
            os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam.bai"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ),
        
        expand(
            os.path.join("COMO_input","{tissue_name}","strandedness","{sample}","{tissue_name}_{tag}_strandedness.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            sample=get.sample(config=config),
            tag=get.tags(config=config),
        ) if perform.get_rnaseq_metrics(config=config) else [],
        
        expand(
            os.path.join(
                config["ROOTDIR"],"data","{tissue_name}","picard","rnaseq","{tissue_name}_{tag}_rnaseq.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ) if perform.get_rnaseq_metrics(config=config) else [],
        
        rule_all_insert_sizes if perform.get_insert_size(config=config) else [],
        expand(os.path.join("COMO_input","{tissue_name}","insertSizeMetrics","{sample}","{tissue_name}_{tag}_insert_size.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            sample=get.sample(config=config)
        ) if perform.get_insert_size(config=config) else [],
        
        rule_all_fragment_sizes if perform.get_fragment_size(config=config) else [],
        expand(
            os.path.join("COMO_input","{tissue_name}","fragmentSizes","{sample}","{tissue_name}_{tag}_fragment_size.txt"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            sample=get.sample(config=config)
        ) if perform.get_fragment_size(config=config) else [],
        
        expand(
            os.path.join(config[
                "ROOTDIR"],"data","{tissue_name}","multiqc",str(config_file_basename),f"{config_file_basename}_multiqc_report.html"),
            tissue_name=get.tissue_name(config=config)
        ),


rule preroundup:
    output:
        layout=os.path.join(config["ROOTDIR"],"data","{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
        preparation=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","prepMethods","{tissue_name}_{tag}_prep_method.txt"),
    resources:
        mem_mb=64,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    run:
        # SRR12873784,effectorcd8_S1R1,PE,total
        sample_row: pd.DataFrame = samples[samples["sample"].str.contains(f"{wildcards.tissue_name}_{wildcards.tag}")]
        name: str = sample_row["sample"].values[0]
        layout: str = sample_row["endtype"].values[0].upper()  # PE, SE, or SLC
        prep_method: str = str(sample_row["prep_method"].values[0].lower())  # total or mrna
        tissue_name: str = name.split("_")[0]
        tag: str = name.split("_")[1]  # S1R1
        study: str = re.match(r"S\d+",tag).group()  # S1
        
        # Write paired/single end or single cell to the appropriate location
        layouts_root: Path = Path(config["ROOTDIR"],"data",tissue_name,"layouts",f"{name}_layout.txt")
        layouts_como: Path = Path("COMO_input",tissue_name,"layouts",study,f"{name}_layout.txt")
        layouts_root.parent.mkdir(parents=True,exist_ok=True)
        layouts_como.parent.mkdir(parents=True,exist_ok=True)
        layouts_write_root = open(layouts_root,"w")
        layouts_write_como = open(layouts_como,"w")
        if Layout[layout] == Layout.PE:
            layouts_write_root.write(Layout.PE.value)
            layouts_write_como.write(Layout.PE.value)
        elif Layout[layout] == Layout.SE:
            layouts_write_root.write(Layout.SE.value)
            layouts_write_como.write(Layout.SE.value)
        elif Layout[layout] == Layout.SLC:
            layouts_write_root.write(Layout.SLC.value)
            layouts_write_como.write(Layout.SLC.value)
        else:
            raise ValueError(f"Invalid selection {layout}. Should be one of 'PE', 'SE', or 'SLC'")
        layouts_write_root.close()
        layouts_write_como.close()
        
        # Write mrna/total to the appropriate location
        prep_root = Path(config["ROOTDIR"],"data",tissue_name,"prepMethods",f"{name}_prep_method.txt")
        prep_como = Path("COMO_input",tissue_name,"prepMethods",study,f"{name}_prep_method.txt")
        prep_root.parent.mkdir(parents=True,exist_ok=True)
        prep_como.parent.mkdir(parents=True,exist_ok=True)
        write_prep_root = open(str(prep_root),"w")
        write_prep_como = open(str(prep_como),"w")
        if PrepMethod[prep_method] == PrepMethod.total:
            write_prep_root.write(PrepMethod.total.value)
            write_prep_como.write(PrepMethod.total.value)
        elif PrepMethod[prep_method] == PrepMethod.mrna or PrepMethod[prep_method] == PrepMethod.polya:
            write_prep_root.write(PrepMethod.mrna.value)
            write_prep_como.write(PrepMethod.mrna.value)
        else:
            raise ValueError(f"Invalid selection {prep_method}. Should be one of 'total', 'mrna', or 'polya'")
        write_prep_root.close()
        write_prep_como.close()
        
        # Make the required directories
        directories: list[str] = [
            os.path.join("COMO_input",tissue_name,"geneCounts"),
            os.path.join("COMO_input",tissue_name,"insertSizeMetrics"),
            os.path.join("COMO_input",tissue_name,"layouts"),
            os.path.join("COMO_input",tissue_name,"layouts",study),
            os.path.join("COMO_input",tissue_name,"fragmentSizes"),
            os.path.join("COMO_input",tissue_name,"prepMethods"),
            os.path.join("COMO_input",tissue_name,"prepMethods",study),
            os.path.join(config["ROOTDIR"],"data",tissue_name,"layouts"),
            os.path.join(config["ROOTDIR"],"data",tissue_name,"prepMethods")
        ]
        for i in directories:
            os.makedirs(name=i,exist_ok=True)


rule generate_genome:
    input:
        genome_fasta_file=config["GENERATE_GENOME"]["GENOME_FASTA_FILE"],
        gtf_file=config["GENERATE_GENOME"]["GTF_FILE"]
    output:
        genome_dir=directory(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"]),
        rule_complete=touch(os.path.join(config["GENERATE_GENOME"]["GENOME_SAVE_DIR"],"generate_genome.complete"))
    threads: 10
    resources:
        mem_mb=51200,
        runtime=150,
        tissue_name="",
    conda: "envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99
        """


rule get_contaminant_genomes:
    output:
        Adapters=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Adapters")),
        Arabidopsis=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Arabidopsis")),
        Drosophila=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Drosophila")),
        E_coli=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"E_coli")),
        Human=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Human")),
        Lambda=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Lambda")),
        Mitochondria=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Mitochondria")),
        Mouse=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Mouse")),
        PhiX=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"PhiX")),
        Rat=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Rat")),
        Vectors=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Vectors")),
        Worm=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Worm")),
        Yeast=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"Yeast")),
        rRNA=directory(os.path.join(config["CONTAMINANT_GENOME_PATH"],"rRNA")),
        config=os.path.join(config["CONTAMINANT_GENOME_PATH"],"fastq_screen.conf")
    threads: 4
    conda: "envs/gnu_parallel.yaml"
    params:
        root_output=config["CONTAMINANT_GENOME_PATH"],
        download_urls=[
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/fastq_screen.conf",
            # These three are the largest genomes, start their download first
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Human",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Mouse",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Rat",
            # --------------------
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Adapters",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Arabidopsis",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Drosophila",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/E_coli",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Lambda",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Mitochondria",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/PhiX",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Vectors",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Worm",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/Yeast",
            "http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/rRNA",
        ],
    resources:
        mem_mb=10240,
        runtime=240,
        tissue_name=""
    shell:
        """
        # Cite parallel before it shows a warning
        echo 'will cite' | parallel --citation > /dev/null 2>&1
        
        download_url() {{
            url="$1"
            
            # Split the url by '/' and get the last item
            item_name=$(echo "$url" | rev | cut -d'/' -f1 | rev)
            
            # Append '/' to the url if the item_name is not 'fastq_screen.conf' to download the entire directory
            [[ "$item_name" != "fastq_screen.conf" ]] && url+="/"
            
            # Test if $item_name exists, only download if it doesn't exist
            if [[ -d "{params.root_output}/$item_name" && "$item_name" != "fastq_screen.conf" ]]; then
                find "{params.root_output}/$item_name" -type f -exec touch {{}} \;
            else
                wget --quiet --recursive --no-parent --no-host-directories --cut-dirs=2 --reject="index.html*" -P {params.root_output} "$url"
            fi
            
            echo "Finished $item_name"
        }}
        
        export -f download_url
        parallel -j {threads} download_url ::: {params.download_urls}
        
        # Wait for all downloads to be done
        wait
        
        # Replace "[FastQ_Screen_Genomes_Path]" with the output directory
        cat "{output.config}" | sed 's|\[FastQ_Screen_Genomes_Path\]|{params.root_output}|g' > "{output.config}.tmp"
        mv "{output.config}.tmp" "{output.config}"
        """


rule prefetch:
    input: config["MASTER_CONTROL"]
    output:
        os.path.join(
            config["ROOTDIR"],"temp","prefetch","{tissue_name}","{tissue_name}_{tag}","{tissue_name}_{tag}.sra")
    conda: "envs/SRAtools.yaml"
    threads: 1
    params:
        srr_value=lambda wildcards: samples.at[
            samples["sample"].eq(f"{wildcards.tissue_name}_{wildcards.tag}").idxmax(),
            "srr"
        ],
        scratch_dir=config["SCRATCH_DIR"],
        temp_file=os.path.join(config['SCRATCH_DIR'],"{tissue_name}_{tag}.sra"),
        output_directory=os.path.join(config["ROOTDIR"],"temp","prefetch","{tissue_name}_{tag}")
    resources:
        mem_mb=16384,
        runtime=20,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","prefetch","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        # If the SRA file lock exists, remove it
        rm -f {output}.lock
        
        # Change into the "scratch" directory so temp files do not populate in the working directory
        curr_dir=$(pwd)
        cd {params.scratch_dir}
        
        prefetch --max-size u --progress --resume yes --output-file {params.temp_file} {params.srr_value}
        
        # Change back to the working directory before moving files
        cd $curr_dir
        mv {params.temp_file} {output}
        
        # Move dependencies into the output directory, checking if files exist in config["SCRATCH_DIR"]
        if [ -n "$(find {params.scratch_dir} -prune -empty)" ]; then
            mv {params.scratch_dir}/* {params.output_directory}
        fi
        """


checkpoint fasterq_dump:
    input:
        prefetch=rules.prefetch.output
    output: fastq=os.path.join(config["ROOTDIR"],"data","{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz")
    threads: 10
    conda: "envs/SRAtools.yaml"
    params:
        scratch_dir=config["SCRATCH_DIR"],
        temp_filename=lambda wildcards: f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq" if wildcards.PE_SE in [
            "1", "2"]
        else f"{wildcards.tissue_name}_{wildcards.tag}.fastq",
        gzip_file=lambda wildcards: f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq.gz" if wildcards.PE_SE in [
            "1", "2"]
        else f"{wildcards.tissue_name}_{wildcards.tag}.fastq.gz",
        split_files=lambda wildcards: True if wildcards.PE_SE in ["1", "2"] else False
    resources:
        mem_mb=10240,
        runtime=45,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","fasterq_dump","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        command='fasterq-dump --force --progress --threads {threads} --temp {params.scratch_dir} --outdir {params.scratch_dir}'
        
        # Set the split/concatenate based on paired end or single end data
        if [[ "{params.split_files}" == "True" ]]; then
            command+=' --split-files'
        else
            command+=' --concatenate-reads'
        fi
        
        # Add the SRA file path to the command
        command+=' {input.prefetch}'
        
        echo $command
        eval $command
        
        # gzip the output
        pigz --synchronous --processes {threads} --force {params.scratch_dir}/{params.temp_filename}
        
        mv {params.scratch_dir}/{params.gzip_file} {output}
        """


def fastqc_dump_fastq_input(wildcards):
    """
    This function will return the input for fastqc_dump_fastq
    It is going to return forward read AND reverse read if the input file is the forward read
    If input is the reverse read, it will only return the reverse read
    If input is a single end read, it will only return the single end read
    """
    if perform.prefetch(config=config):
        if str(wildcards.PE_SE) == "1":
            return [
                checkpoints.fasterq_dump.get(
                    tissue_name=wildcards.tissue_name,
                    tag=wildcards.tag,
                    PE_SE="1"
                ).output[0],
                checkpoints.fasterq_dump.get(
                    tissue_name=wildcards.tissue_name,
                    tag=wildcards.tag,
                    PE_SE="2"
                ).output[0]
            ]
        return checkpoints.fasterq_dump.get(**wildcards).output
    
    for path, subdir, files in os.walk(config["LOCAL_FASTQ_FILES"]):
        for file in files:
            if (
                (wildcards.tissue_name in file) and
                (wildcards.tag in file) and
                (f"_{wildcards.PE_SE}" in file)
            ):
                file_one: str = str(os.path.join(path,file))
                return [file_one,
                        file_one.replace("_1.fastq.gz","_2.fastq.gz")
                        ] if str(wildcards.PE_SE) == "1" else file_one


rule fastqc_dump_fastq:
    input:
        fastq=fastqc_dump_fastq_input,
    output:
        os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params:
        file_one_zip=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        file_one_html=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.html"),
        file_two_zip=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_2_fastqc.zip"),
        file_two_html=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_2_fastqc.html"),
        
        file_one_zip_rename=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"),
        file_one_html_rename=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.html"),
        file_two_zip_rename=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_2_fastqc.zip"),
        file_two_html_rename=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_2_fastqc.html")
    threads: 8
    conda: "envs/fastqc.yaml"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","fastqc_dump_fastq","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
            mv "{params.file_two_zip}" "{params.file_two_zip_rename}"
            mv "{params.file_two_html}" "{params.file_two_html_rename}"

        # Touch the output file because Snakemake will complain about missing files otherwise
        # This file will be created when PE_SE == "1"
        elif [ "{wildcards.PE_SE}" == "2" ]; then
            touch {output}

        elif [ "{wildcards.PE_SE}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{params.file_one_zip_rename}"
            mv "{params.file_one_html}" "{params.file_one_html_rename}"
        fi
        """


def contaminant_screen_input(wildcards):
    # If not performing contaminant screening, return empty list
    if not perform.screen(config=config):
        return []
    
    if perform.trim(config=config):
        return checkpoints.trim.get(**wildcards).output
    
    # If we have performed fasterq_dump, return its output
    if perform.dump_fastq(config=config):
        return checkpoints.fasterq_dump.get(**wildcards).output
    
    # Otherwise collect local files
    fastq_files = Path(config["LOCAL_FASTQ_FILES"])
    for file in fastq_files.rglob("{tissue_name}_{tag}_{PE_SE}.fastq.gz".format(**wildcards)):
        return str(file)


def get_contaminant_genomes_directories(wildcards):
    if not perform.screen(config=config):
        return []
    return contaminant_screen_directories


rule contaminant_screen:
    input:
        files=contaminant_screen_input,
        genomes=get_contaminant_genomes_directories
    output:
        os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","fq_screen","{tissue_name}_{tag}_{PE_SE}_screen.txt")
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        PE_SE="{PE_SE}",
        genomes_config=rules.get_contaminant_genomes.output.config,
        output_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}","fq_screen")
    conda: "envs/screen.yaml"
    threads: 10
    resources:
        mem_mb=6144,
        runtime=30,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","contaminant_screen","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {params.genomes_config} --outdir {params.output_directory} {input.files}
        """


def get_trim_input(wildcards):
    if perform.dump_fastq(config=config):
        if str(wildcards.PE_SE) in ["1", "2"]:
            return [
                checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output,
                checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output
            ]
        return checkpoints.fasterq_dump.get(**wildcards).output
    else:
        pattern = "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
        file = str(next(Path(config["LOCAL_FASTQ_FILES"]).rglob(pattern.format(**wildcards))))
        if str(wildcards.PE_SE) == "1":
            return file, file.replace("_1.fastq.gz","_2.fastq.gz")
        elif str(wildcards.PE_SE) == "2":
            return file.replace("_2.fastq.gz","_1.fastq.gz"), file
        else:
            return file


checkpoint trim:
    input: get_trim_input,
    output:
        os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz")
    threads: 16
    conda: "envs/trim.yaml"
    params:
        scratch_dir=config["SCRATCH_DIR"],
    resources:
        mem_mb=4096,
        runtime=15,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","trim","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        # Find the path to cutadapt, otherwise trim_galore (and it's underlying cutadapt usage) will not work
        cutadapt_path=$(find .snakemake/conda -type f -name "cutadapt")
        
        if [[ "{wildcards.PE_SE}" == "1" ]]; then
            trim_galore --paired --cores 4 -o {params.scratch_dir} --path_to_cutadapt "$cutadapt_path" {input}
            mv "{params.scratch_dir}/{wildcards.tissue_name}_{wildcards.tag}_1_val_1.fq.gz" "{output}"
        elif [[ "{wildcards.PE_SE}" == "2" ]]; then
            trim_galore --paired --cores 4 -o {params.scratch_dir} --path_to_cutadapt "$cutadapt_path" {input}
            mv "{params.scratch_dir}/{wildcards.tissue_name}_{wildcards.tag}_2_val_2.fq.gz" "{output}"
        elif [[ "{wildcards.PE_SE}" == "S" ]]; then
            trim_galore --cores 4 -o {params.scratch_dir} --path_to_cutadapt "$cutadapt_path" {input}
            mv "{params.scratch_dir}/{wildcards.tissue_name}_{wildcards.tag}_S_trimmed.fq.gz" "{output}"
        fi
        """


def get_fastqc_trim_input(wildcards):
    if wildcards.PE_SE == "1":
        forward: str = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output)
        reverse: str = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output)
        return [forward, reverse]
    else:
        return checkpoints.trim.get(**wildcards).output


rule fastqc_trim:
    input: get_fastqc_trim_input  # Original: rules.trim.output
    output:
        os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip")
    params:
        file_two_input=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_2.fastq.gz"),
        file_two_out=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_2_fastqc.zip")
    threads: 8
    conda: "envs/fastqc.yaml"
    resources:
        mem_mb=10240,
        runtime=150,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","fastqc_trim","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

        # Skip reverse reads, but create the output file so Snakemake does not complain about missing files
        # This file will be created when wildcards.PE_SE == "1"
        elif [ "{wildcards.PE_SE}" == "2" ]; then
            touch {output}

        elif [ "{wildcards.PE_SE}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"
        fi
        """


def collect_star_align_input(wildcards):
    if perform.trim(config=config):
        # Have not expanded output from rule trim, need to expand it here
        in_files = sorted(
            expand(
                rules.trim.output,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config)
            )
        )
    elif perform.dump_fastq(config=config):
        # already expanding output from dump_fastq, no need to expand it here
        in_files = sorted(
            expand(
                rules.fasterq_dump.output,
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config)
            )
        )
    
    grouped_reads = []
    for i, in_file in enumerate(in_files):
        direction = get.direction_from_name(in_file)
        try:
            next_file = in_files[i + 1]
            next_direction = get.direction_from_name(next_file)
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
            grouped_reads.append(in_file)
        elif direction == "1" and next_direction == "2":
            if in_file[:-10] == next_file[:-10]:  # remove _1.fastq.gz to make sure they are same replicate
                both_reads = " ".join([in_file, next_file])
                grouped_reads.append(both_reads)
            else:
                warnings.warn(f"{in_file} and {next_file} are incorrectly called together, either the file order is getting scrambled or one end of {in_file} and one end of {next_file} failed to download")
        
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
    This is much like what was done in the function get_dump_fastq_output, located above rule all
    """
    for read in grouped_reads:
        if wildcards.tissue_name in read and wildcards.tag in read:
            return read.split(" ")


def new_star_input(wildcards):
    sample_name: str = f"{wildcards.tissue_name}_{wildcards.tag}"
    is_paired_end: bool = "PE" in samples[samples["sample"].str.startswith(sample_name)]["endtype"].tolist()
    file_pattern = "_1.fastq.gz" if is_paired_end else "_S.fastq.gz"
    pe_suffix = "1" if is_paired_end else "S"
    
    if perform.trim(config=config):
        return [checkpoints.trim.get(**wildcards,PE_SE=pe_suffix).output,
                checkpoints.trim.get(**wildcards,PE_SE="2").output]
    elif perform.dump_fastq(config=config):
        return [checkpoints.fasterq_dump.get(**wildcards,PE_SE=pe_suffix).output,
                checkpoints.fasterq_dump.get(**wildcards,PE_SE="2").output, ]
    else:
        for file in Path(config["LOCAL_FASTQ_FILES"]).rglob(f"*{file_pattern}"):
            if wildcards.tissue_name in str(file) and wildcards.tag in str(file):
                if is_paired_end:
                    print(str(file),str(file).replace(file_pattern,"_2.fastq.gz"))
                    return str(file), str(file).replace(file_pattern,"_2.fastq.gz")
                return str(file)


rule star_align:
    input:
        # reads=collect_star_align_input,
        reads=new_star_input,
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
        gene_table_output=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "aligned_reads",
            "{tag}",
            "{tissue_name}_{tag}_ReadsPerGene.out.tab"
        ),
        bam_output=os.path.join(
            config["ROOTDIR"],
            "data",
            "{tissue_name}",
            "aligned_reads",
            "{tag}",
            "{tissue_name}_{tag}_Aligned.sortedByCoord.out.bam"
        ),
        prefix=os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_"),
    threads: 15
    conda: "envs/star.yaml"
    resources:
        mem_mb=40960,
        runtime=60,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","star_align","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        STAR \
        --runThreadN {threads} \
		--readFilesCommand "zcat" \
		--readFilesIn {input.reads} \
		--genomeDir {input.genome_dir} \
		--outFileNamePrefix {params.prefix} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--quantMode GeneCounts
		
		mv {params.gene_table_output} {output.gene_table}
		mv {params.bam_output} {output.bam_file}
        """

rule index_bam_file:
    input: rules.star_align.output.bam_file
    output: os.path.join(config["ROOTDIR"],"data","{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam.bai")
    threads: 10
    resources:
        mem_mb=1024,
        runtime=10,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda: "envs/samtools.yaml"
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","index_bam_file","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
       samtools index -@ {threads} {input} {output}
       """


rule get_rnaseq_metrics:
    input:
        bam=rules.star_align.output.bam_file,
        tab=rules.star_align.output.gene_table
    output:
        metrics=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","picard","rnaseq","{tissue_name}_{tag}_rnaseq.txt"),
        strand=os.path.join(config["ROOTDIR"],"data","{tissue_name}","strand","{tissue_name}_{tag}_strand.txt")
    params:
        ref_flat=config["GENERATE_GENOME"]["REF_FLAT_FILE"],
        ribo_int_list=config["GENERATE_GENOME"]["RRNA_INTERVAL_LIST"]
    threads: 4
    resources:
        mem_mb=6144,
        runtime=90,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda: "envs/picard.yaml"
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","get_rnaseq_metrics","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        # Create the parent output directories
        #mkdir -p $(dirname -- "{output.metrics}")
        #mkdir -p $(dirname -- "{output.strand}")
        
        # Get the column sums and store them in unst, forw, and rev, respectively
        # We are interested in columns 2, 3, and 4, which correspond to the number of reads in the unstranded, forward, and reverse strand, respectively
        # Column 1: Gene ID
        # Column 2: Counts for unstranded RNA-seq
        # Column 3: Counts for  1st read strand aligned with RNA
        # Column 4:Counts for 2nd read strand aligned with RNA
        colsums=$(grep -v "N_" {input.tab} | awk '{{unstranded+=$2;forward+=$3;reverse+=$4}}END{{print unstranded,forward,reverse}}') || colsums="0 1 2"
        
        # Split colsums based on space (create array of three items)
        IFS=" "
        read -ra arr <<< "$colsums"
        
        # Declare these variables as integers
        declare -i unstranded=${{arr[0]}}
        declare -i forward=${{arr[1]}}
        declare -i reverse=${{arr[2]}}
        
        # Increment the denominator by 1 to prevent "divide by 0"
        if [[ $(( reverse / (forward+1) )) -gt 2 ]]; then
            strand_spec="SECOND_READ_TRANSCRIPTION_STRAND"
        elif [[ $(( forward / (reverse+1) )) -gt 2 ]]; then
            strand_spec="FIRST_READ_TRANSCRIPTION_STRAND"
        else
            strand_spec="NONE"
        fi
        
        echo $strand_spec > {output.strand}
        
        picard CollectRnaSeqMetrics \
            I={input.bam} \
            O={output.metrics} \
            REF_FLAT={config[GENERATE_GENOME][REF_FLAT_FILE]} \
            STRAND_SPECIFICITY=$strand_spec \
            RIBOSOMAL_INTERVALS={config[GENERATE_GENOME][RRNA_INTERVAL_LIST]}
        """


rule get_insert_size:
    input:
        bam=rules.star_align.output.bam_file,
        preround=rules.preroundup.output.layout
    output:
        txt=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","picard","insert","{tissue_name}_{tag}_insert_size.txt"),
        pdf=os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","picard","hist","{tissue_name}_{tag}_insert_size_histo.pdf"),
    threads: 4
    resources:
        mem_mb=1024,
        runtime=60,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda: "envs/picard.yaml"
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","get_insert_size","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        lay=$(cat {input.preround})
        if [ $lay == "paired-end"]; then
            picard CollectinsertSizeMetrics \
            I={input.bam} \
            O={output.txt} \
            H={output.pdf} \
            M=0.05
        else
            echo "cannot collect metrics for single-end data" > {output.txt}
            touch {output.pdf}
        fi
        
        
        
        """

rule get_fragment_size:
    input:
        bam=rules.star_align.output.bam_file,
        bai=rules.index_bam_file.output
    output:
        os.path.join(config["ROOTDIR"],"data","{tissue_name}","fragmentSizes","{tissue_name}_{tag}_fragment_length.txt")
    params:
        layout=os.path.join(config["ROOTDIR"],"data","{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
    threads: 1
    resources:
        partition="batch",
        mem_mb=1024,
        runtime=120,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda: "envs/rseqc.yaml"
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","get_fragment_size","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        # get matches of script file
        file_path=$(find .snakemake/conda/*/bin/RNA_fragment_size.py)
        python3 $file_path -r {config[GENERATE_GENOME][BED_FILE]} -i {input.bam} > {output}
        """

rule copy_gene_counts:
    input: rules.star_align.output.gene_table
    output: os.path.join("COMO_input","{tissue_name}","geneCounts","{sample}","{tissue_name}_{tag}.tab")
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell: """cp {input} {output}"""


rule copy_rnaseq_metrics:
    input: rules.get_rnaseq_metrics.output.strand
    output: os.path.join("COMO_input","{tissue_name}","strandedness","{sample}","{tissue_name}_{tag}_strandedness.txt")
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell: """cp {input} {output}"""

rule copy_insert_size:
    input: rules.get_insert_size.output.txt
    output: os.path.join("COMO_input","{tissue_name}","insertSizeMetrics","{sample}","{tissue_name}_{tag}_insert_size.txt")
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell: """cp {input} {output}"""

rule copy_fragment_size:
    input: rules.get_fragment_size.output
    output: os.path.join("COMO_input","{tissue_name}","fragmentSizes","{sample}","{tissue_name}_{tag}_fragment_size.txt")
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell: """cp {input} {output}"""


def multiqc_get_dump_fastq_data(wildcards):
    if perform.prefetch(config=config):
        output = expand(
            os.path.join(
                config["ROOTDIR"],
                "data",
                "{tissue_name}",
                "raw",
                "{tissue_name}_{tag}_{PE_SE}.fastq.gz"
            ),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        return list(Path(config["LOCAL_FASTQ_FILES"]).rglob(f"*{wildcards.tissue_name}*"))


def multiqc_get_fastqc_data(wildcards):
    if perform.trim(config=config):
        output_files = expand(
            rules.fastqc_trim.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
        )
    else:
        output_files = expand(
            rules.fastqc_dump_fastq.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
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
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


def multiqc_get_screen_data(wildcards):
    if perform.screen(config=config):
        output_files = expand(
            rules.contaminant_screen.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
            PE_SE=get.PE_SE(config=config)
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
    if perform.get_insert_size(config=config):
        for file in expand(
            rules.get_insert_size.output.txt,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_fragmentsize_data(wildcards):
    return_files = []
    if perform.get_fragment_size(config=config):
        for file in expand(
            rules.get_fragment_size.output,
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config)
        ):
            if wildcards.tissue_name in file:
                return_files.append(file)
    return return_files


def multiqc_get_rnaseq_data(wildcards):
    return_files = []
    for file in expand(
        rules.get_rnaseq_metrics.output,
        zip,
        tissue_name=get.tissue_name(config=config),
        tag=get.tags(config=config)
    ):
        if wildcards.tissue_name in file:
            return_files.append(file)
    return return_files


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
        output_file=os.path.join(config[
            "ROOTDIR"],"data","{tissue_name}","multiqc",str(config_file_basename),f"{config_file_basename}_multiqc_report.html"),
        output_directory=directory(os.path.join(
            config["ROOTDIR"],"data","{tissue_name}","multiqc",str(config_file_basename)))
    params:
        config_file_basename=config_file_basename,
        input_directory=os.path.join(config["ROOTDIR"],"data","{tissue_name}")
    threads: 1
    conda: "envs/multiqc.yaml"
    resources:
        mem_mb=5120,
        runtime=30,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(os.path.join("benchmarks","{tissue_name}","multiqc","{tissue_name}.benchmark"),
            config["BENCHMARK_TIMES"])
    shell:
        """
        mkdir -p "{output.output_directory}"
        
        multiqc --interactive --force \
            --title "{wildcards.tissue_name}" \
            -filename {params.config_file_basename}_multiqc_report.html \
            --outdir "{output.output_directory}" \
            "{params.input_directory}"
        
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
