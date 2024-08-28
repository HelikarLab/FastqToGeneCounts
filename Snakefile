import csv
import os
import pandas as pd
import warnings
from pathlib import Path

from utils import get, perform
from utils.constants import Layout, PrepMethod
from utils.genome_generation import Utilities


configfile: "config.yaml"

os.makedirs(config["ROOTDIR"],exist_ok=True)

# Get the delimiter from the master control file; from: https://stackoverflow.com/questions/16312104
delimiter = csv.Sniffer().sniff(open(config["MASTER_CONTROL"]).read(1024)).delimiter
samples: pd.DataFrame = pd.read_csv(
    filepath_or_buffer=str(config["MASTER_CONTROL"]),
    delimiter=str(delimiter),
    names=["srr", "sample", "endtype", "prep_method"],
)

config_file_basename: str = os.path.basename(config["MASTER_CONTROL"]).split(".")[0]
screen_genomes: pd.DataFrame = pd.read_csv("utils/screen_genomes.csv",delimiter=",",header=0)
contaminant_genomes_root = os.path.join(config["ROOTDIR"],"FastQ_Screen_Genomes")
root_data = os.path.join(config["ROOTDIR"],"data")
root_temp = os.path.join(config["ROOTDIR"],"temp")
species_name = Utilities.get_species_from_taxon(taxon_id=config["GENOME"]["TAXONOMY_ID"])

if config["GENOME"]["VERSION"] == "latest":
    ensembl_release_number = f"release-{Utilities.get_latest_release()}"
elif config["GENOME"]["VERSION"].startswith("release"):
    ensembl_release_number = config["GENOME"]["VERSION"]
elif config["GENOME"]["VERSION"].isdigit():
    ensembl_release_number = f"release-{config['GENOME']['VERSION']}"
else:
    raise ValueError(
        "Invalid GENOME VERSION in config.yaml file. Valid options are: 'latest', 'release-###' (i.e., 'release-112'), or an integer (i.e., 112)"
    )


rule all:
    input:
        # Genome generation items + star genome index
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}.bed"),
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_genome_sizes.txt"),
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_{ensembl_release_number}.gtf"),
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_rrna.interval_list"),
        os.path.join(
            config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_{ensembl_release_number}_primary_assembly.fa"
        ),
        os.path.join(
            config["GENOME"]["SAVE_DIR"],
            species_name,
            f"{species_name}_{ensembl_release_number}_primary_assembly.fa.fai",
        ),
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_ref_flat.txt"),
        os.path.join(config["GENOME"]["SAVE_DIR"],species_name,"star","job_complete.txt"),
        expand(
            os.path.join(root_data,"{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        ),
        expand(
            os.path.join(root_data,"{tissue_name}","prepMethods","{tissue_name}_{tag}_prep_method.txt"),
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        ),
        expand(
            os.path.join(
                root_data,
                "{tissue_name}",
                "fastqc",
                "untrimmed_reads",
                "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
            ),
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            PE_SE=get.PE_SE(config),
        ),
        expand(
            os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab"),
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        ),
        expand(
            os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam"),
            zip,
            tissue_name=get.tissue_name(config=config),
            tag=get.tags(config=config),
        ),
        expand(
            os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam.bai"),
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        ),
        expand(
            os.path.join(
                root_data,
                "{tissue_name}",
                "multiqc",
                config_file_basename,
                f"{config_file_basename}_multiqc_report.html",
            ),
            tissue_name=get.tissue_name(config),
        ),
        (
            expand(
                os.path.join(root_temp,"prefetch","{tissue_name}","{tissue_name}_{tag}","{tissue_name}_{tag}.sra"),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
            )
            if perform.prefetch(config=config)
            else []
        ),
        (
            expand(
                os.path.join(config["ROOTDIR"],"FastQ_Screen_Genomes","{name}"),
                name=screen_genomes["name"].values,
            )
            if perform.screen(config=config)
            else []
        ),
        (
            expand(
                os.path.join(root_data,"{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config),
            )
            if perform.prefetch(config=config) or perform.dump_fastq(config=config)
            else []
        ),
        (
            expand(
                os.path.join(
                    root_data,"{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"
                ),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config),
            )
            if perform.trim(config=config)
            else []
        ),
        (
            expand(
                os.path.join(root_data,"{tissue_name}","fragmentSizes","{tissue_name}_{tag}_fragment_length.txt"),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
            )
            if perform.get_fragment_size(config=config)
            else []
        ),
        (
            expand(
                os.path.join(root_data,"{tissue_name}","fq_screen","{tissue_name}_{tag}_{PE_SE}_screen.txt"),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config),
            )
            if perform.screen(config=config)
            else []
        ),
        (
            expand(
                os.path.join(
                    root_data,
                    "{tissue_name}",
                    "fastqc",
                    "trimmed_reads",
                    "trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip",
                ),
                zip,
                tissue_name=get.tissue_name(config=config),
                tag=get.tags(config=config),
                PE_SE=get.PE_SE(config=config),
            )
            if perform.trim(config=config)
            else []
        ),
        (
            expand(
                os.path.join(root_data,"{tissue_name}","picard","rnaseq","{tissue_name}_{tag}_rnaseq.txt"),
                zip,
                tissue_name=get.tissue_name(config),
                tag=get.tags(config),
            )
            if perform.get_rnaseq_metrics(config)
            else []
        ),
        expand(
            os.path.join("COMO_input","{tissue_name}","geneCounts","{sample}","{tissue_name}_{tag}.tab"),
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            sample=get.sample(config),
        ),
        (
            expand(
                os.path.join(
                    "COMO_input","{tissue_name}","strandedness","{sample}","{tissue_name}_{tag}_strandedness.txt"
                ),
                zip,
                tissue_name=get.tissue_name(config),
                sample=get.sample(config),
                tag=get.tags(config),
            )
            if perform.get_rnaseq_metrics(config)
            else []
        ),
        (
            expand(
                os.path.join(
                    "COMO_input",
                    "{tissue_name}",
                    "insertSizeMetrics",
                    "{sample}",
                    "{tissue_name}_{tag}_insert_size.txt",
                ),
                zip,
                tissue_name=get.tissue_name(config),
                tag=get.tags(config),
                sample=get.sample(config),
            )
            if perform.get_insert_size(config)
            else []
        ),
        (
            expand(
                os.path.join(
                    "COMO_input","{tissue_name}","fragmentSizes","{sample}","{tissue_name}_{tag}_fragment_size.txt"
                ),
                zip,
                tissue_name=get.tissue_name(config),
                tag=get.tags(config),
                sample=get.sample(config),
            )
            if perform.get_fragment_size(config)
            else []
        ),


rule preroundup:
    output:
        layout=os.path.join(root_data,"{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
        preparation=os.path.join(root_data,"{tissue_name}","prepMethods","{tissue_name}_{tag}_prep_method.txt"),
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    run:
        # SRR12873784,effectorcd8_S1R1,PE,total
        sample_row: pd.DataFrame = samples.loc[
                                   samples["sample"] == f"{wildcards.tissue_name}_{wildcards.tag}",
                                   :,
                                   # Collect everything from the row by using `:`
                                   ]
        # Collect the required data
        srr_code: str = sample_row["srr"].values[0]
        name: str = sample_row["sample"].values[0]
        endtype: str = sample_row["endtype"].values[0].upper()
        prep_method: str = sample_row["prep_method"].values[0].lower()
        tissue_name: str = name.split("_")[0]
        tag: str = name.split("_")[1]
        study: str = re.match(r"S\d+",tag).group()

        # Write paired/single end or single cell to the appropriate location
        layouts_root: Path = Path(root_data,tissue_name,"layouts",f"{name}_layout.txt")
        layouts_como: Path = Path("COMO_input",tissue_name,"layouts",study,f"{name}_layout.txt")
        layouts_root.parent.mkdir(parents=True,exist_ok=True)
        layouts_como.parent.mkdir(parents=True,exist_ok=True)
        layouts_write_root = open(layouts_root,"w")
        layouts_write_como = open(layouts_como,"w")
        layout: str = str(sample_row["endtype"].values[0]).upper()  # PE, SE, or SLC

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
        prep_root: Path = Path(root_data,tissue_name,"prepMethods",f"{name}_prep_method.txt")
        prep_como: Path = Path("COMO_input",tissue_name,"prepMethods",study,f"{name}_prep_method.txt")
        prep_root.parent.mkdir(parents=True,exist_ok=True)
        prep_como.parent.mkdir(parents=True,exist_ok=True)
        write_prep_root = open(prep_root.as_posix(),"w")
        write_prep_como = open(prep_como.as_posix(),"w")
        prep_method = str(sample_row["prep_method"].values[0]).lower()  # total or mrna

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
            os.path.join(root_data,tissue_name,"layouts"),
            os.path.join(root_data,tissue_name,"prepMethods"),
            os.path.join("COMO_input",tissue_name,"geneCounts"),
            os.path.join("COMO_input",tissue_name,"insertSizeMetrics"),
            os.path.join("COMO_input",tissue_name,"layouts"),
            os.path.join("COMO_input",tissue_name,"layouts",study),
            os.path.join("COMO_input",tissue_name,"fragmentSizes"),
            os.path.join("COMO_input",tissue_name,"prepMethods"),
            os.path.join("COMO_input",tissue_name,"prepMethods",study),
            os.path.join(root_data,tissue_name,"layouts"),
            os.path.join(root_data,tissue_name,"prepMethods"),
        ]
        for i in directories:
            os.makedirs(name=i,exist_ok=True)


rule generate_genome:
    output:
        bed_file=os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}.bed"),
        ref_flat=os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_ref_flat.txt"),
        genome_sizes=os.path.join(config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_genome_sizes.txt"),
        gtf_file=os.path.join(
            config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_{ensembl_release_number}.gtf"
        ),
        rrna_interval_list=os.path.join(
            config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_rrna.interval_list"
        ),
        primary_assembly=os.path.join(
            config["GENOME"]["SAVE_DIR"],species_name,f"{species_name}_{ensembl_release_number}_primary_assembly.fa"
        ),
        primary_assembly_index=os.path.join(
            config["GENOME"]["SAVE_DIR"],
            species_name,
            f"{species_name}_{ensembl_release_number}_primary_assembly.fa.fai",
        ),
    conda:
        "envs/generate_genome.yaml"
    threads: 1
    resources:
        mem_mb=8096,
        runtime=30,
        tissue_name="",# intentionally left blank, reference: github.com/jdblischak/smk-simple-slurm/issues/20
    shell:
        """
        python3 utils/genome_generation.py \
            --taxon-id {config[GENOME][TAXONOMY_ID]} \
            --release-number {config[GENOME][VERSION]} \
            --root-save-dir {config[GENOME][SAVE_DIR]}
        """


rule star_index_genome:
    input:
        primary_assembly=rules.generate_genome.output.primary_assembly,
        gtf_file=rules.generate_genome.output.gtf_file,
    output:
        job_complete=os.path.join(config["GENOME"]["SAVE_DIR"],species_name,"star","job_complete.txt"),
    params:
        species_name=species_name,
    conda:
        "envs/star.yaml"
    threads: 10
    resources:
        mem_mb=51200,
        runtime=150,
        tissue_name="",# intentionally left blank, reference: github.com/jdblischak/smk-simple-slurm/issues/20
    shell:
        """
        genome_directory=$(dirname {output.job_complete})
        mkdir -p "$genome_directory"

        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir "$genome_directory" \
        --genomeFastaFiles {input.primary_assembly} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99

        echo "STAR genomeGenerate for species '{params.species_name}' finished at $(date)" > {output.job_complete}
        """


rule get_contaminant_genomes:
    output:
        root_output=directory(contaminant_genomes_root),
        Adapters=directory(os.path.join(contaminant_genomes_root,"Adapters")),
        Arabidopsis=directory(os.path.join(contaminant_genomes_root,"Arabidopsis")),
        Drosophila=directory(os.path.join(contaminant_genomes_root,"Drosophila")),
        E_coli=directory(os.path.join(contaminant_genomes_root,"E_coli")),
        Human=directory(os.path.join(contaminant_genomes_root,"Human")),
        Lambda=directory(os.path.join(contaminant_genomes_root,"Lambda")),
        Mitochondria=directory(os.path.join(contaminant_genomes_root,"Mitochondria")),
        Mouse=directory(os.path.join(contaminant_genomes_root,"Mouse")),
        PhiX=directory(os.path.join(contaminant_genomes_root,"PhiX")),
        Rat=directory(os.path.join(contaminant_genomes_root,"Rat")),
        Vectors=directory(os.path.join(contaminant_genomes_root,"Vectors")),
        Worm=directory(os.path.join(contaminant_genomes_root,"Worm")),
        Yeast=directory(os.path.join(contaminant_genomes_root,"Yeast")),
        rRNA=directory(os.path.join(contaminant_genomes_root,"rRNA")),
        config=os.path.join(contaminant_genomes_root,"fastq_screen.conf"),
    threads: 1
    params:
        zip_url=r"https://uofnelincoln-my.sharepoint.com/:u:/g/personal/jloecker3_unl_edu/EWO6p5t-kjZEks3pW-1RVvgBR2Q-nFI_8kXx5_NB8_kYnw?e=dp2NWX&download=1",
    resources:
        mem_mb=6144,
        runtime=30,
        tissue_name="",# intentionally left blank, reference: github.com/jdblischak/smk-simple-slurm/issues/20
    shell:
        """
        mkdir -p {output.root_output}
        wget --quiet "{params.zip_url}" -O "{output.root_output}/FastQ_Screen_Genomes.zip"

        # Unzip the archive
        unzip -o "{output.root_output}/FastQ_Screen_Genomes.zip" -d "{output.root_output}"
        rm "{output.root_output}/FastQ_Screen_Genomes.zip"

        # Replace "[FastQ_Screen_Genomes_Path]" with the output directory, then remove any double slashes (//)
        sed 's|\\[FastQ_Screen_Genomes_Path\\]|{output.root_output}|g' "{output.config}" | sed 's|//|/|g' > "{output.config}.tmp"

        mv "{output.config}.tmp" "{output.config}"
        """


rule prefetch:
    input:
        config["MASTER_CONTROL"],
    output:
        os.path.join(root_temp,"prefetch","{tissue_name}","{tissue_name}_{tag}","{tissue_name}_{tag}.sra"),
    conda:
        "envs/SRAtools.yaml"
    threads: 1
    params:
        srr_value=lambda wildcards: samples.at[
            samples["sample"].eq(f"{wildcards.tissue_name}_{wildcards.tag}").idxmax(),
            "srr",
        ],
        scratch_dir=config["SCRATCH_DIR"],
        output_directory=os.path.join(root_temp,"prefetch","{tissue_name}_{tag}"),
        temp_file=os.path.join(config["SCRATCH_DIR"],"{tissue_name}_{tag}.sra"),
    resources:
        mem_mb=16384,
        runtime=20,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","prefetch","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        echo Starting rule
        rm -f {output}.lock

        echo Making scratch directory
        mkdir -p {params.scratch_dir}

        echo Starting prefetch
        prefetch --max-size u --progress --output-file {params.temp_file} {params.srr_value}

        echo Making output directory
        mkdir -p "$(dirname {output})"

        echo Moving files
        mv {params.scratch_dir}/* "$(dirname {output})/"
        """


checkpoint fasterq_dump:
    input:
        prefetch=rules.prefetch.output,
    output:
        fastq=os.path.join(root_data,"{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
    threads: 10
    conda:
        "envs/SRAtools.yaml"
    params:
        scratch_dir=config["SCRATCH_DIR"],
        temp_filename=lambda wildcards: (
            f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq"
            if wildcards.PE_SE in ["1", "2"]
            else f"{wildcards.tissue_name}_{wildcards.tag}.fastq"
        ),
        gzip_file=lambda wildcards: (
            f"{wildcards.tissue_name}_{wildcards.tag}_{wildcards.PE_SE}.fastq.gz"
            if wildcards.PE_SE in ["1", "2"]
            else f"{wildcards.tissue_name}_{wildcards.tag}.fastq.gz"
        ),
        split_files=lambda wildcards: True if wildcards.PE_SE in ["1", "2"] else False,
    resources:
        mem_mb=10240,
        runtime=45,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","fasterq_dump","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        command='fasterq-dump --force --progress --threads {threads} --temp {params.scratch_dir} --outdir {params.scratch_dir}'

        # Set the split/concatenate based on paired end or single end data
        [[ "{params.split_files}" == "True" ]] && command+=' --split-files' || command+=' --concatenate-reads'

        # Add the SRA file path to the command
        command+=' {input.prefetch}'

        eval $command
        ls -l {params.scratch_dir}
        pigz --synchronous --processes {threads} --force {params.scratch_dir}/{params.temp_filename}
        mv {params.scratch_dir}/{params.gzip_file} {output.fastq}
        """


def fastqc_dump_fastq_input(wildcards):
    """
    This function will return the input for fastqc_dump_fastq
    It is going to return forward read AND reverse read if the input file is the forward read
    If input is the reverse read, it will only return the reverse read
    If input is a single end read, it will only return the single end read
    """
    if perform.prefetch(config):
        if str(wildcards.PE_SE) == "1":
            return [
                checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output[0],
                checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output[0],
            ]
        return checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="S").output

    # Make sure we are able to load local FastQ files
    if "LOCAL_FASTQ_FILES" in config.keys() and os.path.exists(config["LOCAL_FASTQ_FILES"]):
        for path, subdir, files in os.walk(config["LOCAL_FASTQ_FILES"]):
            for file in files:
                if (wildcards.tissue_name in file) and (wildcards.tag in file) and (f"_{wildcards.PE_SE}" in file):
                    file_one: str = str(os.path.join(path,file))
                    return (
                        [file_one, file_one.replace("_1.fastq.gz","_2.fastq.gz")]
                        if str(wildcards.PE_SE) == "1"
                        else file_one
                    )


rule fastqc_dump_fastq:
    input:
        fastq=fastqc_dump_fastq_input,
    output:
        file_one_zip_rename=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"
        ),
        file_one_html_rename=os.path.join(
            root_data,
            "{tissue_name}",
            "fastqc",
            "untrimmed_reads",
            "untrimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.html",
        ),
    params:
        file_one_zip=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.zip"
        ),
        file_one_html=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_{PE_SE}_fastqc.html"
        ),
        file_two_zip=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_2_fastqc.zip"
        ),
        file_two_html=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","{tissue_name}_{tag}_2_fastqc.html"
        ),
        file_two_zip_rename=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_2_fastqc.zip"
        ),
        file_two_html_rename=os.path.join(
            root_data,"{tissue_name}","fastqc","untrimmed_reads","untrimmed_{tissue_name}_{tag}_2_fastqc.html"
        ),
    threads: 8
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","fastqc_dump_fastq","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        output_directory="$(dirname {output.file_one_zip_rename})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{output.file_one_zip_rename}"
            mv "{params.file_one_html}" "{output.file_one_html_rename}"
            mv "{params.file_two_zip}" "{params.file_two_zip_rename}"
            mv "{params.file_two_html}" "{params.file_two_html_rename}"

        # Touch the output file because Snakemake will complain about missing files otherwise
        # This file will be created when PE_SE == "1"
        elif [ "{wildcards.PE_SE}" == "2" ]; then
            touch {params.file_two_zip_rename}
            touch {params.file_two_html_rename}

        elif [ "{wildcards.PE_SE}" == "S" ]; then
            fastqc {input} --threads {threads} -o "$output_directory"

            mv "{params.file_one_zip}" "{output.file_one_zip_rename}"
            mv "{params.file_one_html}" "{output.file_one_html_rename}"
        fi
        """


def contaminant_screen_input(wildcards):
    # If not performing contaminant screening, return empty list
    if not perform.screen(config):
        return []

    # If we have performed fasterq_dump, return its output
    if perform.dump_fastq(config):
        return checkpoints.fasterq_dump.get(**wildcards).output

    # Otherwise collect local files
    fastq_files = Path(config["LOCAL_FASTQ_FILES"])
    for file in fastq_files.rglob("{tissue_name}_{tag}_{PE_SE}.fastq.gz".format(**wildcards)):
        return str(file)


rule contaminant_screen:
    input:
        files=contaminant_screen_input,
        genomes=contaminant_genomes_root,
    output:
        os.path.join(root_data,"{tissue_name}","fq_screen","{tissue_name}_{tag}_{PE_SE}_screen.txt"),
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        PE_SE="{PE_SE}",
        genomes_config=rules.get_contaminant_genomes.output.config,
        output_directory=os.path.join(root_data,"{tissue_name}","fq_screen"),
    conda:
        "envs/screen.yaml"
    threads: 10
    resources:
        mem_mb=6144,
        runtime=30,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","contaminant_screen","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {params.genomes_config} --outdir {params.output_directory} {input.files}
        """


def get_trim_input(wildcards):
    if perform.dump_fastq(config):
        if str(wildcards.PE_SE) in ["1", "2"]:
            return [
                *checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output,
                *checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output,
            ]
        return checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="S").output
    else:
        files = list(
            Path(config["LOCAL_FASTQ_FILES"]).rglob("{tissue_name}_{tag}_{PE_SE}.fastq.gz".format(**wildcards))
        )
        if str(wildcards.PE_SE) in ["1", "2"]:
            return [str(file) for file in files] + [str(file).replace("_1.fastq.gz","_2.fastq.gz") for file in files]
        else:
            return [str(file) for file in files]


checkpoint trim:
    input:
        get_trim_input,
    output:
        os.path.join(root_data,"{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
    threads: 16
    conda:
        "envs/trim.yaml"
    params:
        scratch_dir=config["SCRATCH_DIR"],
    resources:
        mem_mb=10240,
        runtime=120,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","trim","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        output_directory="$(dirname {output})"

        if [[ "{wildcards.PE_SE}" == "1" ]]; then
            file_out_1="{params.scratch_dir}/{wildcards.tissue_name}_{wildcards.tag}_1_val_1.fq.gz"    # final output paired end, forward read
            trim_galore --paired --cores 4 -o {params.scratch_dir} {input}
            mv "$file_out_1" "{output}"

        # Skip over reverse-reads. Create the output file so snakemake does not complain about the rule not generating output
        elif [[ "{wildcards.PE_SE}" == "2" ]]; then
            file_out_2="{params.scratch_dir}/{wildcards.tissue_name}_{wildcards.tag}_2_val_2.fq.gz"    # final output paired end, reverse read
            trim_galore --paired --cores 4 -o {params.scratch_dir} {input}
            mv "$file_out_2" "{output}"

        # Work on single-end reads
        elif [[ "{wildcards.PE_SE}" == "S" ]]; then
            file_out="$output_directory/{wildcards.tissue_name}_{wildcards.tag}_S_trimmed.fq.gz"   # final output single end
            trim_galore --cores 4 -o "$output_directory" {input}
            mv "$file_out" "{output}"
        fi
        """


def get_fastqc_trim_input(wildcards):
    if wildcards.PE_SE == "1":
        forward = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output)
        reverse = str(checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output)
        return [forward, reverse]
    else:
        return checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="S").output


rule fastqc_trim:
    input:
        get_fastqc_trim_input,
    output:
        os.path.join(
            root_data,"{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_{PE_SE}_fastqc.zip"
        ),
    params:
        file_two_input=os.path.join(
            root_data,"{tissue_name}","trimmed_reads","trimmed_{tissue_name}_{tag}_2.fastq.gz"
        ),
        file_two_out=os.path.join(
            root_data,"{tissue_name}","fastqc","trimmed_reads","trimmed_{tissue_name}_{tag}_2_fastqc.zip"
        ),
    threads: 8
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=10240,
        runtime=150,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","fastqc_trim","{tissue_name}_{tag}_{PE_SE}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        output_directory="$(dirname {output})"
        mkdir -p "$output_directory"

        if [ "{wildcards.PE_SE}" == "1" ]; then
            # send fastqc commands to background so we can run both at the same time
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
    rule_output = rules.trim.output if perform.trim else rules.fasterq_dump.output
    in_files = sorted(
        expand(
            rule_output,
            zip,
            tissue_name=wildcards.tissue_name,
            tag=wildcards.tag,
            PE_SE=get.PE_SE(config),
        )
    )

    grouped_reads = []
    for i, in_file in enumerate(in_files):
        direction = direction_from_name(in_file)
        try:
            next_file = in_files[i + 1]
            next_direction = direction_from_name(next_file)
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
    This is much like what was done in the function get_dump_fastq_output, located above rule all
    """
    for read in grouped_reads:
        if wildcards.tissue_name in read and wildcards.tag in read:
            return read.split(" ")


def new_star_input(wildcards):
    items = []
    sample_name: str = f"{wildcards.tissue_name}_{wildcards.tag}"
    is_paired_end: bool = "PE" == samples[samples["sample"] == sample_name]["endtype"].values[0]
    file_pattern = "_1.fastq.gz" if is_paired_end else "_S.fastq.gz"
    if perform.trim(config):
        items.extend(
            [
                *checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output,
                *checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output,
            ]
            if is_paired_end
            else [*checkpoints.trim.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="S").output]
        )
    elif perform.dump_fastq(config):
        items.extend(
            [
                *checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="1").output,
                *checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="2").output,
            ]
            if is_paired_end
            else [*checkpoints.fasterq_dump.get(tissue_name=wildcards.tissue_name,tag=wildcards.tag,PE_SE="S").output]
        )
    else:
        for file in Path(config["LOCAL_FASTQ_FILES"]).rglob(f"*{file_pattern}"):
            if wildcards.tissue_name in file.as_posix() and wildcards.tag in file.as_posix():
                items.extend(
                    [file.as_posix(), file.as_posix().replace(file_pattern,"_2.fastq.gz")]
                    if is_paired_end
                    else file.as_posix()
                )
                break

    return items


rule star_align:
    input:
        reads=new_star_input,
        genome_index=rules.star_index_genome.output.job_complete,
    output:
        gene_table=os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.tab"),
        bam_file=os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam"),
    params:
        tissue_name="{tissue_name}",
        tag="{tag}",
        gene_table_output=os.path.join(
            root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_ReadsPerGene.out.tab"
        ),
        bam_output=os.path.join(
            root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_Aligned.sortedByCoord.out.bam"
        ),
        prefix=os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}_"),
    threads: 15
    conda:
        "envs/star.yaml"
    resources:
        mem_mb=40960,
        runtime=60,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","star_align","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        genome_directory=$(dirname {input.genome_index})

        STAR \
        --runThreadN {threads} \
        --readFilesCommand "zcat" \
        --readFilesIn {input.reads} \
        --genomeDir "$genome_directory" \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts

        mv {params.gene_table_output} {output.gene_table}
        mv {params.bam_output} {output.bam_file}
        """


rule index_bam_file:
    input:
        rules.star_align.output.bam_file,
    output:
        os.path.join(root_data,"{tissue_name}","aligned_reads","{tag}","{tissue_name}_{tag}.bam.bai"),
    conda:
        "envs/samtools.yaml"
    threads: 10
    resources:
        mem_mb=1024,
        runtime=10,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","index_bam_file","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
       samtools index -@ {threads} {input} {output}
       """


rule get_rnaseq_metrics:
    input:
        bam=rules.star_align.output.bam_file,
        tab=rules.star_align.output.gene_table,
        ref_flat=rules.generate_genome.output.ref_flat,
        rrna_interval_list=rules.generate_genome.output.rrna_interval_list,
    output:
        metrics=os.path.join(root_data,"{tissue_name}","picard","rnaseq","{tissue_name}_{tag}_rnaseq.txt"),
        strand=os.path.join(root_data,"{tissue_name}","strand","{tissue_name}_{tag}_strand.txt"),
    threads: 4
    resources:
        mem_mb=6144,
        runtime=90,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda:
        "envs/picard.yaml"
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","get_rnaseq_metrics","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        # Create the parent output directories
        #mkdir -p $(dirname -- "{output.metrics}")

        # Get the column sums and store them in unst, forw, and rev, respectively
        # We are interested in columns 2, 3, and 4, which correspond to the number of reads in the unstranded, forward, and reverse strand, respectively
        # Column 1: Gene ID
        # Column 2: Counts for unstranded RNA-seq
        # Column 3: Counts for  1st read strand aligned with RNA
        # Column 4: Counts for 2nd read strand aligned with RNA
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

        picard CollectRnaSeqMetrics I={input.bam} O={output.metrics} REF_FLAT={input.ref_flat} STRAND_SPECIFICITY=$strand_spec RIBOSOMAL_INTERVALS={input.rrna_interval_list}
        """


rule get_insert_size:
    input:
        bam=rules.star_align.output.bam_file,
        preround=rules.preroundup.output.layout,
    output:
        txt=os.path.join(root_data,"{tissue_name}","picard","insert","{tissue_name}_{tag}_insert_size.txt"),
        pdf=os.path.join(root_data,"{tissue_name}","picard","hist","{tissue_name}_{tag}_insert_size_histo.pdf"),
    threads: 4
    resources:
        mem_mb=1024,
        runtime=60,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda:
        "envs/picard.yaml"
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","get_insert_size","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        lay=$(cat {input.preround})
        if [ $lay == "paired-end"]; then
            picard CollectinsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf} M=0.05
            picard CollectinsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf} M=0.05
        else
            echo "cannot collect metrics for single-end data" > {output.txt}
            touch {output.pdf}
        fi
        """


rule get_fragment_size:
    input:
        bam=rules.star_align.output.bam_file,
        bai=rules.index_bam_file.output,
        bed_file=rules.generate_genome.output.bed_file,
    output:
        os.path.join(root_data,"{tissue_name}","fragmentSizes","{tissue_name}_{tag}_fragment_length.txt"),
    params:
        layout=os.path.join(root_data,"{tissue_name}","layouts","{tissue_name}_{tag}_layout.txt"),
    threads: 1
    resources:
        partition="batch",
        mem_mb=1024,
        runtime=120,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    conda:
        "envs/rseqc.yaml"
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","get_fragment_size","{tissue_name}_{tag}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        # get matches of script file
        file_path=$(find .snakemake/conda/*/bin/RNA_fragment_size.py)
        python3 $file_path -r {input.bed_file} -i {input.bam} > {output}
        """


rule copy_gene_counts:
    input:
        rules.star_align.output.gene_table,
    output:
        os.path.join("COMO_input","{tissue_name}","geneCounts","{sample}","{tissue_name}_{tag}.tab"),
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell:
        """cp --link {input} {output}"""


rule copy_rnaseq_metrics:
    input:
        rules.get_rnaseq_metrics.output.strand,
    output:
        os.path.join("COMO_input","{tissue_name}","strandedness","{sample}","{tissue_name}_{tag}_strandedness.txt"),
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell:
        """cp --link {input} {output}"""


rule copy_insert_size:
    input:
        rules.get_insert_size.output.txt,
    output:
        os.path.join(
            "COMO_input","{tissue_name}","insertSizeMetrics","{sample}","{tissue_name}_{tag}_insert_size.txt"
        ),
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell:
        """cp --link {input} {output}"""


rule copy_fragment_size:
    input:
        rules.get_fragment_size.output,
    output:
        os.path.join(
            "COMO_input","{tissue_name}","fragmentSizes","{sample}","{tissue_name}_{tag}_fragment_size.txt"
        ),
    threads: 1
    resources:
        mem_mb=256,
        runtime=1,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    shell:
        """cp --link {input} {output}"""


def multiqc_get_dump_fastq_data(wildcards):
    return (
        expand(
            os.path.join(root_data,"{tissue_name}","raw","{tissue_name}_{tag}_{PE_SE}.fastq.gz"),
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            PE_SE=get.PE_SE(config),
        )
        if perform.prefetch(config)
        else list(Path(config["LOCAL_FASTQ_FILES"]).rglob(f"*{wildcards.tissue_name}*"))
    )


def multiqc_get_fastqc_data(wildcards):
    output_files = (
        expand(
            rules.fastqc_trim.output,
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            PE_SE=get.PE_SE(config),
        )
        if perform.trim(config)
        else expand(
            rules.fastqc_dump_fastq.output,
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            PE_SE=get.PE_SE(config),
        )
    )
    return_files = [file for file in output_files if wildcards.tissue_name in file]
    return return_files


def multiqc_get_star_data(wildcards):
    output_files = expand(
        rules.star_align.output.gene_table,
        zip,
        tissue_name=get.tissue_name(config),
        tag=get.tags(config),
    )
    return_files = [file for file in output_files if wildcards.tissue_name in file]
    return return_files


def multiqc_get_screen_data(wildcards):
    output_files = (
        expand(
            rules.contaminant_screen.output,
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
            PE_SE=get.PE_SE(config),
        )
        if perform.screen(config)
        else []
    )
    return_files = [file for file in output_files if wildcards.tissue_name in file]
    return return_files


def multiqc_get_insertsize_data(wildcards):
    output_files = (
        expand(
            rules.get_insert_size.output.txt,
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        )
        if perform.get_insert_size(config)
        else []
    )
    return_files = (
        [file for file in output_files if wildcards.tissue_name in file] if perform.get_insert_size(config) else []
    )
    return return_files


def multiqc_get_fragmentsize_data(wildcards):
    output_files = (
        expand(
            rules.get_fragment_size.output,
            zip,
            tissue_name=get.tissue_name(config),
            tag=get.tags(config),
        )
        if perform.get_fragment_size(config)
        else []
    )
    return_files = (
        [file for file in output_files if wildcards.tissue_name in file] if perform.get_fragment_size(config) else []
    )
    return return_files


def multiqc_get_rnaseq_data(wildcards):
    output_files = expand(
        rules.get_rnaseq_metrics.output,
        zip,
        tissue_name=get.tissue_name(config),
        tag=get.tags(config),
    )
    return_files = [file for file in output_files if wildcards.tissue_name in file]
    return return_files


rule multiqc:
    input:
        fastqc_data=multiqc_get_fastqc_data,
        star_data=multiqc_get_star_data,
        dump_fastq_data=multiqc_get_dump_fastq_data,
        screen_data=multiqc_get_screen_data,
        insertsize_data=multiqc_get_insertsize_data,
        rnaseq_data=multiqc_get_rnaseq_data,
        fragment_size_data=multiqc_get_fragmentsize_data,
    output:
        output_file=os.path.join(
            root_data,
            "{tissue_name}",
            "multiqc",
            str(config_file_basename),
            f"{config_file_basename}_multiqc_report.html",
        ),
        output_directory=directory(os.path.join(root_data,"{tissue_name}","multiqc",str(config_file_basename))),
    params:
        config_file_basename=config_file_basename,
        input_directory=os.path.join(root_data,"{tissue_name}"),
    threads: 1
    conda:
        "envs/multiqc.yaml"
    resources:
        mem_mb=5120,
        runtime=30,
        tissue_name=lambda wildcards: wildcards.tissue_name,
    benchmark:
        repeat(
            os.path.join("benchmarks","{tissue_name}","multiqc","{tissue_name}.benchmark"),
            config["BENCHMARK_TIMES"],
        )
    shell:
        """
        mkdir -p "{output.output_directory}"

        multiqc --interactive --force --title "{wildcards.tissue_name}" --filename {params.config_file_basename}_multiqc_report.html --outdir "{output.output_directory}" "{params.input_directory}"
        """
