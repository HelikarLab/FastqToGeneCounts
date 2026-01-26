import os
import re
import sys
from pathlib import Path
from typing import Literal

import pandas as pd

from utils.constants import Layout, PrepMethod
from utils.parse import Config, SampleData, print_key_value_table

configfile: "config.yaml"
cfg: Config = Config.create(config)
data: SampleData = SampleData(cfg.sample_filepath)

# Only print the table if we are not calculating the dag (for png/svg creation, etc.) at the start
onstart:
    if "--rulegraph" not in sys.argv and "--dag" not in sys.argv:
        print_key_value_table(
            "Parsed Filepaths & Directories",
            [
                ("Samples", cfg.sample_filepath),
                ("Data", cfg.data_root),
                ("Temporary", cfg.temp_root),
                ("COMO", cfg.como_root),
                ("Logging", cfg.logs_root),
                ("Genome", cfg.genome.species_dir),
            ],
        )

rule all:
    input:
        f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_primary_assembly.fa",
        f"{cfg.genome.contaminants_dir}/.download_complete",

        f"{cfg.genome.species_dir}/star/Log.out",
        f"{cfg.genome.species_dir}/transcriptome.fa",

        expand(f"{cfg.data_root}/{{tissue}}/layouts/{{tissue}}_{{tag}}_layout.txt", zip, tissue=data.tissues, tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/prepMethods/{{tissue}}_{{tag}}_prep_method.txt",zip,tissue=data.tissues,tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam",zip,tissue=data.tissues,tag=data.tags),
        expand(f"{cfg.como_root}/{{tissue}}/geneCounts/{{study}}/{{tissue}}_{{tag}}.tab",zip,tissue=data.tissues,study=data.studies,tag=data.tags),
        expand(f"{cfg.como_root}/{{tissue}}/strandedness/{{study}}/{{tissue}}_{{tag}}_strandedness.txt",zip,tissue=data.tissues,study=data.studies,tag=data.tags),
        expand(f"{cfg.como_root}/{{tissue}}/insertSizeMetrics/{{study}}/{{tissue}}_{{tag}}_insert_size.txt",zip,tissue=data.tissues,study=data.studies,tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}/{cfg.sample_filepath.stem}_multiqc_report.html",tissue=data.tissues),
        expand(
            f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_{{end}}_fastqc.zip",
            zip,
            tissue=data.tissues_paired,
            tag=data.tags_paired,
            end=data.ends_paired
        ),
        expand(
            f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_{{end}}_fastqc.zip",
            zip,
            tissue=data.tissues_paired,
            tag=data.tags_paired,
            end=data.ends_paired
        ),

        branch(
            perform.dump_fastq(config),
            then=expand(f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
            otherwise=[],
        ),
        branch(
            perform.trim(config),
            then=expand(f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
            otherwise=[],
        ),
        branch(
            perform.screen(config),
            then=f"{cfg.genome.contaminants_dir}/fastq_screen.conf",
            otherwise=[],
        ),
        branch(
            perform.get_fragment_size(config),
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/fragmentSizes/{{tissue}}_{{tag}}_fragment_size.txt", zip, tissue=data.tissues, tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/fragmentSizes/{{study}}/{{tissue}}_{{tag}}_fragment_size.txt",zip,tissue=data.tissues,study=data.studies,tag=data.tags)
            ],
            otherwise=[]
        ),
        branch(
            perform.get_rnaseq_metrics(config),
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/picard/rnaseq/{{tissue}}_{{tag}}_rnaseq.txt", zip, tissue=data.tissues, tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/strandedness/{{study}}/{{tissue}}_{{tag}}_strandedness.txt",zip,tissue=data.tissues,tag=data.tags,study=data.studies),
            ],
            otherwise=[],
        ),
        branch(
            perform.get_insert_size(config),
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.txt",zip,tissue=data.tissues,tag=data.tags),
                expand(f"{cfg.data_root}/{{tissue}}/picard/hist/{{tissue}}_{{tag}}_insert_size_histo.pdf",zip,tissue=data.tissues,tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/insertSizeMetrics/{{study}}/{{tissue}}_{{tag}}_insert_size.txt",zip,tissue=data.tissues,tag=data.tags,study=data.studies),
            ],
            otherwise=[]
        )

rule preroundup:
    output:
        layout=f"{cfg.data_root}/{{tissue}}/layouts/{{tissue}}_{{tag}}_layout.txt",
        preparation=f"{cfg.data_root}/{{tissue}}/prepMethods/{{tissue}}_{{tag}}_prep_method.txt",
    params:
        sample_name=lambda wildcards: f"{wildcards.tissue}_{wildcards.tag}",
    resources:
        mem_mb=1024,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    run:
        # example row: SRR12873784,effectorcd8_S1R1,PE,total
        sample_row: pd.Series = samples[samples["sample"].eq(params.sample_name)]
        endtype: str = sample_row["endtype"].values[0].upper()
        prep_method: str = sample_row["prep_method"].values[0].lower()
        study = re.match(r"S\d+",wildcards.tag).group()

        # Make the required directories
        (cfg.data_root / wildcards.tissue / "layouts").mkdir(parents=True,exist_ok=True)
        (cfg.data_root / wildcards.tissue / "prepMethods").mkdir(parents=True,exist_ok=True)
        (cfg.como_root / wildcards.tissue / "layouts" / study).mkdir(parents=True,exist_ok=True)
        (cfg.como_root / wildcards.tissue / "prepMethods" / study).mkdir(parents=True,exist_ok=True)

        # Write paired/single end or single cell to the appropriate location
        layouts_root: Path = Path(cfg.data_root,wildcards.tissue,"layouts",f"{params.sample_name}_layout.txt")
        layouts_como: Path = Path(cfg.como_root,wildcards.tissue,"layouts",study,f"{params.sample_name}_layout.txt")
        with layouts_root.open("w") as layouts_write_root, layouts_como.open("w") as layouts_write_como:
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

        # Write mrna/total to the appropriate location
        prep_root: Path = Path(cfg.data_root,wildcards.tissue,"prepMethods",f"{params.sample_name}_prep_method.txt")
        prep_como: Path = Path(cfg.como_root,wildcards.tissue,"prepMethods",study,f"{params.sample_name}_prep_method.txt")
        with prep_root.open("w") as write_prep_root, prep_como.open("w") as write_prep_como:
            prep_method = str(sample_row["prep_method"].values[0]).lower()  # total or tissue
            if PrepMethod[prep_method] == PrepMethod.total:
                write_prep_root.write(PrepMethod.total.value)
                write_prep_como.write(PrepMethod.total.value)
            elif PrepMethod[prep_method] == PrepMethod.mrna or PrepMethod[prep_method] == PrepMethod.polya:
                write_prep_root.write(PrepMethod.mrna.value)
                write_prep_como.write(PrepMethod.mrna.value)
            else:
                raise ValueError(f"Invalid selection {prep_method}. Should be one of 'total', 'mrna', or 'polya'")

rule download_genome:
    output:
        bed_file=f"{cfg.genome.species_dir}/{cfg.species_name}.bed",
        ref_flat=f"{cfg.genome.species_dir}/{cfg.species_name}_ref_flat.txt",
        genome_sizes=f"{cfg.genome.species_dir}/{cfg.species_name}_genome_sizes.txt",
        gtf_file=f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}.gtf",
        rrna_interval_list=f"{cfg.genome.species_dir}/{cfg.species_name}_rrna.interval_list",
        primary_assembly=f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_primary_assembly.fa",
        primary_assembly_index=f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_primary_assembly.fa.fai",
    conda: "envs/generate_genome.yaml"
    threads: 1
    resources:
        mem_mb=8096,
        runtime=30,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
        network_slots=1
    shell:
        """
        python3 utils/download_genome.py \
            --taxon-id {cfg.genome.taxon_id} \
            --release-number {cfg.genome.version} \
            --root-save-dir {cfg.genome.species_dir}
        """

rule download_contaminant_genomes:
    output:
        done=f"{cfg.genome.contaminants_dir}/.download_complete",
        config=f"{cfg.genome.contaminants_dir}/fastq_screen.conf",
    conda:
        "envs/screen.yaml"
    threads: 1
    params:
        root_output=directory(cfg.genome.contaminants_dir),
    resources:
        mem_mb=6144,
        runtime=30,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
        network_slots=1
    log: f"{cfg.logs_root}/download_contaminant_genomes.log"
    shell:
        r"""
        echo "" > {log}
        
        # get the line that does not have Bisulfite in it
        genome_location=$(curl --silent "https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt" | grep -v "Bisulfite" | head -n 1)
        if [[ "$genome_location" != http*://* ]]; then
            genome_location="https://$genome_location"
        fi

        # Normalize genome location by removing (optional) trailing slash
        genome_base="${{genome_location%/}}"
        mkdir -p "{params.root_output}"

        expected_genomes=("Adapters" "Arabidopsis" "Drosophila" "E_coli" "Human" "Lambda" "Mitochondria" "Mouse" "PhiX" "Rat" "Vectors" "Worm" "Yeast" "rRNA")
        expected_files=(6 6 6 6 6 6 6 6 6 6 6 6 6 7)
        for index in ${{!expected_genomes[@]}}; do
            genome="${{expected_genomes[index]}}"
            outdir="{params.root_output}/$genome"

            # skip genome if directory exists and is non-empty
            if [[ -d "$outdir" ]]; then
                existing_files=$(ls -1 $outdir | wc -l)
                expected_files="${{expected_files[index]}}"
                if [ $existing_files -eq $expected_files ]; then
                    echo "[fastq_screen] Skipping genome download for '$genome' because it is already present at '$outdir'" >> {log}
                    continue
                fi
            fi

            echo "[fastq_screen] Downloading genome '$genome' to: $outdir" >> {log}
            mkdir -p "$outdir"
            wget --quiet --recursive --no-parent --no-host-directories --cut-dirs=4 --reject "index.html*,robots.txt*" --directory-prefix "{params.root_output}" "$genome_base/$genome/"
        done

        if [[ ! -f "{output.config}" ]]; then
            wget --quiet -O "{output.config}" "${{genome_base}}/fastq_screen.conf"
        fi

        # Replace "[FastQ_Screen_Genomes_Path]" with the output directory, then remove any double slashes (//)
        sed "s|\[FastQ_Screen_Genomes_Path\]|{params.root_output}|g" "{output.config}" | sed "s|//|/|g" > "{output.config}.tmp"
        mv "{output.config}.tmp" {output.config}

        touch "{output.done}"
        """

rule star_index_genome:
    input:
        primary_assembly=rules.download_genome.output.primary_assembly,
        gtf_file=rules.download_genome.output.gtf_file,
    output:
        chromosome_length=f"{cfg.genome.species_dir}/star/chrLength.txt",
        chromosome_name_length=f"{cfg.genome.species_dir}/star/chrNameLength.txt",
        chromosome_name=f"{cfg.genome.species_dir}/star/chrName.txt",
        chromosome_start=f"{cfg.genome.species_dir}/star/chrStart.txt",
        exon_gene_info=f"{cfg.genome.species_dir}/star/exonGeTrInfo.tab",
        exon_info=f"{cfg.genome.species_dir}/star/exonInfo.tab",
        gene_info=f"{cfg.genome.species_dir}/star/geneInfo.tab",
        genome=f"{cfg.genome.species_dir}/star/Genome",
        genome_parameters=f"{cfg.genome.species_dir}/star/genomeParameters.txt",
        log=f"{cfg.genome.species_dir}/star/Log.out",
        suffix_array=f"{cfg.genome.species_dir}/star/SA",
        genome_index=f"{cfg.genome.species_dir}/star/SAindex",
        slice_juntion_db_info=f"{cfg.genome.species_dir}/star/sjdbInfo.txt",
        slice_juntion_db_from_gtf=f"{cfg.genome.species_dir}/star/sjdbList.fromGTF.out.tab",
        slice_juntion_db_list=f"{cfg.genome.species_dir}/star/sjdbList.out.tab",
        transcript_info=f"{cfg.genome.species_dir}/star/transcriptInfo.tab",
    params:
        output_dir=f"{cfg.genome.species_dir}/star",
    conda:
        "envs/star.yaml"
    threads: 10
    resources:
        mem_mb=51200,
        runtime=150,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
    log: f"{cfg.logs_root}/star_index_genome.log"
    benchmark:
        repeat(
            # format is '2025-05-16T112835'
            f"{cfg.benchmark_dir}/star_index_genome/star_index_genome_{time.strftime('%Y-%m-%dT%H%M%S',time.localtime())}.benchmark",
            cfg.benchmark_count,
        )
    shell:
        """
        mkdir -p {params.output_dir}

        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {params.output_dir} \
        --genomeFastaFiles {input.primary_assembly} \
        --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang 99 1>{log} 2>&1
        """

rule generate_transcriptome_fasta:
    input:
        genome=rules.download_genome.output.primary_assembly,
        gtf=rules.download_genome.output.gtf_file
    output:
        transcriptome=f"{cfg.genome.species_dir}/transcriptome.fa"
    conda:
        "envs/gffread.yaml"
    log: f"{cfg.logs_root}/generate_transcriptome_fasta.log"
    resources:
        mem_mb=4096,
        time=10,
    shell: "gffread -w {output.transcriptome} -g {input.genome} {input.gtf} 1>{log} 2>&1"

rule prefetch:
    input:
        cfg.sample_filepath,
    output:
        temp(f"{cfg.data_root}/{{tissue}}/.prefetch/{{tissue}}_{{tag}}.sra"),
    conda:
        "envs/SRAtools.yaml"
    params:
        output_directory=f"{cfg.data_root}/prefetch/{{tissue}}_{{tag}}",
        srr=lookup(query="sample == '{tissue}_{tag}'",within=data.samples,cols="srr"),
    threads: 1
    log: f"{cfg.logs_root}/{{tissue}}/prefetch/{{tissue}}_{{tag}}_prefetch.log"
    resources:
        mem_mb=16384,
        runtime=20,
        tissue=lambda wildcards: wildcards.tissue,
        network_slots=1
    shell:
        """
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        sra_temp="$tmpdir/{wildcards.tissue}_{wildcards.tag}.sra"

        prefetch --max-size u --progress --output-file "$sra_temp" {params.srr} 1>{log} 2>&1
        mv "$sra_temp" {output}
        """

rule fastq_dump_paired:
    input:
        rules.prefetch.output
    output:
        r1=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_1.fastq.gz",
        r2=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_2.fastq.gz",
    threads: 5
    log: f"{cfg.logs_root}/{{tissue}}/fastq_dump/{{tissue}}_{{tag}}_fastq_dump.log"
    resources:
        mem_mb=10240,
        runtime=45,
        tissue=lambda wildcards: wildcards.tissue
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        
        tmp_forward="$tmpdir/{wildcards.tissue}_{wildcards.tag}_1.fastq"
        tmp_reverse="$tmpdir/{wildcards.tissue}_{wildcards.tag}_2.fastq"
        fasterq-dump --force --split-files --progress --threads {threads} --temp "$tmpdir" --outdir "$tmpdir" {input} 1>{log} 2>&1
        pigz --processes {threads} --force "$tmp_forward" "$tmp_reverse"
        
        printf "\n\n" >> {log}
        echo "Moving '$tmp_forward' to '{output.r1}'" >> {log}
        echo "Moving '$tmp_reverse' to '{output.r2}'" >> {log}
        mv "$tmp_forward.gz" "{output.r1}" &
        mv "$tmp_reverse.gz" "{output.r2}" &
        
        wait
        """

rule fastq_dump_single:
    input:
        rules.prefetch.output
    output:
        S=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_S.fastq.gz"
    threads: 4
    resources:
        mem_mb=10240,
        runtime=30,
        tissue=lambda wildcards: wildcards.tissue
    log: f"{cfg.logs_root}/{{tissue}}/fastq_dump/{{tissue}}_{{tag}}_fastq_dump.log"
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        
        tmpfile="$tmpdir/{wildcards.tissue}_{wildcards.tag}_S.fastq"
        fasterq-dump --force --concatenate-reads --progress --threads {threads} --temp "$tmpdir" --outdir "$tmpdir" {input} 1>{log} 2>&1
        pigz --processes 4 --force "$tmpfile"
        
        printf "\n\n" >> {log}
        echo "Moving '$tmpfile' to '{output.S}'" >> {log}
        mv "$tmpfile" {output.S}
        """


def qc_raw_fastq_paired_input(wildcards):
    if perform.prefetch(config):
        # We will get data from NCBI, therefore we can get results from fastq_dump_paired directly
        return [rules.fastq_dump_paired.output.r1, rules.fastq_dump_paired.output.r2]

    # Otherwise, we must search the local directory for fiels and return those
    if "LOCAL_FASTQ_FILES" in config and exists(cfg.local_fastq_filepath):
        sample_name = f"{wildcards.tissue}_{wildcards.tag}"
        for path, subdir, files in os.walk(cfg.local_fastq_filepath):
            for file in files:
                if sample_name in data.sample_names and file.startswith("_".join(wildcards)):
                    forward_read = str(path / file)
                    reverse_read = forward_read.replace("_1.fastq.gz", "_2.fastq.gz")
                    return [forward_read, reverse_read]
    else:
        print("Unable to find directory 'LOCAL_FASTQ_FILES' defined in 'config.yaml'.")
    return []

rule qc_raw_fastq_paired:
    input:
        reads=qc_raw_fastq_paired_input
    output:
        r1_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_1_fastqc.zip",
        r1_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_1_fastqc.html",
        r2_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_2_fastqc.zip",
        r2_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_2_fastqc.html",
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/raw/{{tissue}}_{{tag}}_fastqc.log"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue=lambda wildcards: wildcards.tissue
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input.reads} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip' to '{output.r1_zip}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html' to '{output.r1_html}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip' to '{output.r2_zip}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html' to '{output.r2_html}'" >> {log}
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip" "{output.r1_zip}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html" "{output.r1_html}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip" "{output.r2_zip}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html" "{output.r2_html}"
        """

rule qc_raw_fastq_single:
    input:
        rules.fastq_dump_single.output,
    output:
        s_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_S_fastqc.zip",
        s_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_S_fastqc.html"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/raw/{{tissue}}_{{tag}}_fastqc.log"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue=lambda wildcards: wildcards.tissue
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads 5 -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.zip' to '{output.s_zip}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.html' to '{output.s_html}'" >> {log}
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.zip" "{output.s_zip}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.html" "{output.s_html}"
        """

rule trim_paired:
    input:
        r1=rules.fastq_dump_paired.output.r1,
        r2=rules.fastq_dump_paired.output.r2
    output:
        r1=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_1.fastq.gz",
        r2=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_2.fastq.gz",
    # See the trim_galore `--cores` setting for details on why 16 was chosen
    # https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
    threads: 16
    log: f"{cfg.logs_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_trim.log"
    resources:
        mem_mb=10240,
        runtime=120,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        trim_galore --paired --cores 4 -o "$tmpdir" {input} 1>{log} 2>&1
        
        printf "\n\n" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_val_1.fq.gz' to '{output.r1}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_val_2.fq.gz' to '{output.r2}'" >> {log}
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_val_1.fq.gz" "{output.r1}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_val_2.fq.gz" "{output.r2}"
        """

rule trim_single:
    input:
        rules.fastq_dump_single.output.S
    output:
        S=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_S.fastq.gz"
    # See the trim_galore `--cores` setting for details on why 16 was chosen
    # https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
    threads: 16
    log: f"{cfg.logs_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_trim.log"
    resources:
        mem_mb=10240,
        runtime=120,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        trim_galore --cores 4 -o "$tmpdir" {input} 1>{log} 2>&1
        
        printf "\n\n" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_trimmed.fq.gz' to '{output.S}'" >> {log}
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_trimmed.fq.gz" "{output.S}"
        """

rule qc_trim_fastq_paired:
    input:
        r1=rules.trim_paired.output.r1,
        r2=rules.trim_paired.output.r2
    output:
        r1_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_1_fastqc.zip",
        r1_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_1_fastqc.html",
        r2_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_2_fastqc.zip",
        r2_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_2_fastqc.html",
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/trimmed/{{tissue}}_{{tag}}_fastqc.log"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue=lambda wildcards: wildcards.tissue
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip' to '{output.r1_zip}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html' to '{output.r1_html}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip' to '{output.r2_zip}'" >> {log}
        echo "Moving '$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html' to '{output.r2_html}'" >> {log}
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip" "{output.r1_zip}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html" "{output.r1_html}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip" "{output.r2_zip}"
        mv "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html" "{output.r2_html}"
        """

rule qc_trim_fastq_single:
    input:
        rules.trim_single.output,
    output:
        s_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_S_fastqc.zip",
        s_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_S_fastqc.html"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/trimmed/{{tissue}}_{{tag}}_fastqc.log"
    resources:
        mem_mb=4096,
        runtime=150,
        tissue=lambda wildcards: wildcards.tissue
    shell:
        r"""
        tmpdir=$(mktemp -d)
        tmp_zip="$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.zip"
        tmp_html="$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.html"
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        echo "Moving '$tmp_zip' to '{output.s_zip}'" >> {log}
        echo "Moving '$tmp_html' to '{output.s_html}'" >> {log}
        mv "$tmp_zip" "{output.s_zip}"
        mv "$tmp_html" "{output.s_html}"
        """

rule align:
    input:
        files=lambda wildcards: (
            rules.fastq_dump_paired.output
            if all(data.samples.loc[data.samples["sample"] == f"{wildcards.tissue}_{wildcards.tag}", "endtype"] == "PE")
            else rules.fastq_dump_single.output
        ),
        genome=rules.download_genome.output.primary_assembly
    output:  # Not all outputs listed are used, but they are listed so Snakemake knows to clean them up on workflow reruns
        bam_file=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam",
        final_log=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.final.out",
        intermediate_log=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.out",
        progress=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.progress.out",
        splice_junctions=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_SJ.out.tab",
        gene_table=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.tab",
        tx_bam=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.toTranscriptome.bam"
    conda:
        "envs/star.yaml"
    params:
        star_genome=f"{cfg.genome.species_dir.as_posix().removesuffix('/')}/star",
        gene_table=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_ReadsPerGene.out.tab",
        bam_output=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Aligned.sortedByCoord.out.bam",
        tx_bam_output=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Aligned.toTranscriptome.out.bam",
    threads: 5
    log: f"{cfg.logs_root}/{{tissue}}/align/{{tissue}}_{{tag}}_star_align.log"
    resources:
        mem_mb=32768,
        runtime=15,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        r"""
        # remove any files not listed in the output
        rm -rf "$(dirname {output.gene_table})/*"
        
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        
        STAR \
        --runThreadN {threads} \
        --readFilesCommand "zcat" \
        --readFilesIn {input.files} \
        --genomeDir "{params.star_genome}" \
        --outFileNamePrefix "$tmpdir/{wildcards.tissue}_{wildcards.tag}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts TranscriptomeSAM 1>{log} 2>&1
        
        printf "\n\n" >> {log}
        echo "Moving files from temporary directory '$tmpdir' to '$(dirname {output.gene_table})'" >> {log}
        mv $tmpdir/* "$(dirname {output.gene_table})/"
        
        
        echo "Moving '{params.gene_table}' to '{output.gene_table}'" >> {log}
        echo "Moving '{params.bam_output}' to '{output.bam_file}'" >> {log}
        echo "Moving '{params.tx_bam_output}' to '{output.tx_bam}'" >> {log}
        mv {params.gene_table} {output.gene_table}
        mv {params.bam_output} {output.bam_file}
        mv {params.tx_bam_output} {output.tx_bam}
        """

rule index_bam_file:
    input:
        rules.align.output.bam_file,
    output:
        f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam.bai",
    conda:
        "envs/samtools.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/align/{{tissue}}_{{tag}}_index_bam.log"
    resources:
        mem_mb=1024,
        runtime=10,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/index_bam_file/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """
        samtools index -@ {threads} {input} {output} 1>{log} 2>&1
        """

rule salmon_quantification:
    input:
        tx_bam=rules.align.output.tx_bam,
        transcriptome=rules.generate_transcriptome_fasta.output.transcriptome
    output:
        quant=f"{cfg.data_root}/{{tissue}}/read_quantification/salmon/{{tag}}/{{tissue}}_{{tag}}_quant.sf",
        meta=f"{cfg.data_root}/{{tissue}}/read_quantification/salmon/{{tag}}/{{tissue}}_{{tag}}_meta_info.json"
    params:
        outdir=f"{cfg.data_root}/{{tissue}}/read_quantification/salmon/{{tag}}"
    conda:
        "envs/salmon.yaml"
    threads: 8
    log: f"{cfg.logs_root}/{{tissue}}/salmon_quant/{{tissue}}_{{tag}}_salmon_quant.log"
    resources:
        mem_mb=8192,
        runtime=20,
        tissue=lambda wildcards: wildcards.tissue
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/salmon_quant/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        mkdir -p {params.outdir}

        salmon quant \
            --threads {threads} \
            --targets {input.transcriptome} \
            --libType A \
            --alignments {input.tx_bam} \
            --output {params.outdir} \
            --seqBias --gcBias --posBias --useVBOpt 1>{log} 2>&1

        printf "\n\n" >> {log}
        echo "Moving '{params.outdir}/quant.sf' to '{output.quant}'" >> {log}
        echo "Moving '{params.outdir}/cmd_info.json' to '{output.meta}'" >> {log}
        mv {params.outdir}/quant.sf {output.quant}
        mv {params.outdir}/cmd_info.json {output.meta}
        """

def contaminant_screen_input(wildcards):
    if not perform.screen(config):
        return []

    if perform.dump_fastq(config):
        if all(data.samples.loc[data.samples["sample"] == f"{wildcards.tissue}_{wildcards.tag}", "endtype"] == "PE"):
            return rules.fastq_dump_paired.output
        else:
            return rules.fastq_dump_single.output

    files: list[str] = []
    for file in cfg.local_fastq_filepath.rglob("{tissue}_{tag}_{PE_SE}.fastq.gz".format(**wildcards)):
        if f"{wildcards.tissue}_{wildcards.tag}" in data.samples["sample"].tolist():
            files.append(file.as_posix())
    return files

rule contaminant_screen_paired:
    input:
        files=contaminant_screen_input,
        screen_config = rules.download_contaminant_genomes.output.config,
    output:
        r1=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_1_screen.txt",
        r2=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_2_screen.txt"
    params:
        output_dir=f"{cfg.data_root}/{{tissue}}/fq_screen",
    conda:
        "envs/screen.yaml"
    threads: 5
    log: f"{cfg.logs_root}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_fastq_screen.log"
    resources:
        mem_mb=6144,
        runtime=30,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(
            f"{cfg.benchmark_dir}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_paired.benchmark",
            cfg.benchmark_count
        )
    shell:
        """
        outdir=$(dirname {output.r1})
        mkdir -p "$outdir"
        fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {input.screen_config} --outdir "$outdir" {input.files} 1>{log} 2>&1
        """

rule contaminant_screen_single:
    input:
        files=contaminant_screen_input,
        screen_config=rules.download_contaminant_genomes.output.config
    output:
        S=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_S_screen.txt"
    params:
        output_dir=f"{cfg.data_root}/{{tissue}}/fq_screen"
    conda:
        "envs/screen.yaml"
    threads: 5
    log: f"{cfg.logs_root}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_fastq_screen.log"
    resources:
        mem_mb=6144,
        runtime=30,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_single.benchmark", cfg.benchmark_count)
    shell:
        """
        outdir=$(dirname {output.S})
        mkdir -p "$outdir"
        fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {input.screen_config} --outdir "$outdir" {input.files} 1>{log} 2>&1
        """

rule fragment_size:
    input:
        bam = rules.align.output.bam_file,
        bai = rules.index_bam_file.output,
        bed_file = rules.download_genome.output.bed_file,
    output:
        f"{cfg.data_root}/{{tissue}}/fragmentSizes/{{tissue}}_{{tag}}_fragment_size.txt",
    params:
        layout=f"{cfg.data_root}/{{tissue}}/layouts/{{tissue}}_{{tag}}_layout.txt",
        bed_filepath=rules.download_genome.output,
    conda: "envs/rseqc.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fragment_size/{{tissue}}_{{tag}}_fragment_size.log"
    resources:
        mem_mb=1024,
        runtime=120,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/get_fragment_size/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        "python3 utils/get_fragment_size.py --input {input.bam} --bai {input.bai} --refgene {input.bed_file} --log {log} --output {output}"


rule insert_size:
    input:
        bam=rules.align.output.bam_file,
        layout=rules.preroundup.output.layout
    output:
        txt=f"{cfg.data_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.txt",
        pdf=f"{cfg.data_root}/{{tissue}}/picard/hist/{{tissue}}_{{tag}}_insert_size_histo.pdf",
    conda: "envs/picard.yaml"
    threads: 1
    log: f"{cfg.logs_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.log"
    resources:
        mem_mb=1024,
        runtime=5,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/get_insert_size/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """
        layout=$(cat {input.layout})
        if [ $layout == "paired-end" ]; then
            picard CollectInsertSizeMetrics \
                --VERBOSITY WARNING \
                --INPUT {input.bam} \
                --OUTPUT {output.txt} \
                --Histogram_FILE {output.pdf} \
                --MINIMUM_PCT 0.05 1>{log} 2>&1
        else
            echo "NOTICE: Cannot collect metrics for single-end data." | tee -a "{log}" "{output.txt}"
            touch {output.pdf}
        fi
        """

rule rnaseq_metrics:
    input:
        bam = rules.align.output.bam_file,
        tab = rules.align.output.gene_table,
        ref_flat = rules.download_genome.output.ref_flat,
        rrna_interval_list = rules.download_genome.output.rrna_interval_list,
    output:
        metrics=f"{cfg.data_root}/{{tissue}}/picard/rnaseq/{{tissue}}_{{tag}}_rnaseq.txt",
        strand=f"{cfg.data_root}/{{tissue}}/strand/{{tissue}}_{{tag}}_strand.txt",
    conda: "envs/picard.yaml"
    threads: 1
    log: f"{cfg.logs_root}/{{tissue}}/picard/rnaseq/{{tissue}}_{{tag}}_rnaseq_metrics.log"
    resources:
        mem_mb=1024,
        runtime=5,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/get_rnaseq_metrics/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """
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

        picard CollectRnaSeqMetrics \
            --VERBOSITY WARNING \
            --INPUT {input.bam} \
            --OUTPUT {output.metrics} \
            --REF_FLAT {input.ref_flat} \
            --STRAND_SPECIFICITY $strand_spec \
            --RIBOSOMAL_INTERVALS {input.rrna_interval_list} 1>{log} 2>&1
        """

rule copy_fragment_size:
    localrule: True
    input:
        rules.fragment_size.output,
    output:
        f"{cfg.como_root}/{{tissue}}/fragmentSizes/{{sample}}/{{tissue}}_{{tag}}_fragment_size.txt"
    resources:
        mem_mb=512,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        """cp {input} {output}"""

rule copy_insert_size:
    localrule: True
    input:
        rules.insert_size.output.txt,
    output:
        f"{cfg.como_root}/{{tissue}}/insertSizeMetrics/{{sample}}/{{tissue}}_{{tag}}_insert_size.txt"
    resources:
        mem_mb=512,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        """cp {input} {output}"""

rule copy_rnaseq_metrics:
    localrule: True
    input:
        rules.rnaseq_metrics.output.strand,
    output:
        f"{cfg.como_root}/{{tissue}}/strandedness/{{sample}}/{{tissue}}_{{tag}}_strandedness.txt"
    resources:
        mem_mb=512,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        """cp {input} {output}"""


rule copy_gene_counts:
    localrule: True
    input:
        rules.align.output.gene_table,
    output:
        f"{cfg.como_root}/{{tissue}}/geneCounts/{{sample}}/{{tissue}}_{{tag}}.tab"
    resources:
        mem_mb=512,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    shell:
        """cp {input} {output}"""


rule multiqc:
    input:
        raw_fastq=expand(f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
        trimmed_fastq=expand(f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
        aligned_fastq=expand(rules.align.output.bam_file, zip, tissue=data.tissues, tag=data.tags),
        contamination_paired=expand(rules.contaminant_screen_paired.output,zip,tissue=data.tissues_paired,tag=data.tags_paired),
        contamination_single=expand(rules.contaminant_screen_single.output,zip,tissue=data.tissues_paired,tag=data.tags_paired),
        insert_sizes=expand(rules.insert_size.output.txt, zip, tissue=data.tissues, tag=data.tags),
        rnaseq_metrics=expand(rules.rnaseq_metrics.output.metrics, zip, tissue=data.tissues, tag=data.tags),
        fragment_sizes=expand(rules.fragment_size.output, zip, tissue=data.tissues, tag=data.tags),
        salmon_quant=expand(rules.salmon_quantification.output.quant, zip, tissue=data.tissues, tag=data.tags),
    output:
        output_file=f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}/{cfg.sample_filepath.stem}_multiqc_report.html",
    params:
        config_file_basename=cfg.sample_filepath.stem,
        input_directory=f"{cfg.data_root}/{{tissue}}",
        output_directory=directory(f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}"),
    conda:
        "envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb=5120,
        runtime=30,
        tissue=lambda wildcards: wildcards.tissue,
    log: f"{cfg.logs_root}/{{tissue}}/multiqc/{{tissue}}_multiqc.log"
    benchmark:
        repeat(f"{cfg.benchmark_dir}/{{tissue}}/multiqc/{{tissue}}.benchmark",cfg.benchmark_count)
    shell:
        """
        filename="$(basename {output})"
        mkdir -p "{params.output_directory}"

        multiqc --interactive --force --title "{wildcards.tissue}" \
          --filename {params.config_file_basename}_multiqc_report.html \
          --outdir "{params.output_directory}" {params.input_directory} 1>{log} 2>&1
        """
