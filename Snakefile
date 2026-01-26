import sys
from typing import Literal

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
        expand(f"{cfg.data_root}/{{tissue}}/{{tissue}}_config.yaml", tissue=set(data.tissues)),
        f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_{cfg.genome.type}.fa",
        f"{cfg.genome.contaminants_dir}/.download_complete",
        f"{cfg.genome.species_dir}/star/Log.out",
        f"{cfg.genome.species_dir}/transcriptome.fa",

        expand(f"{cfg.data_root}/{{tissue}}/layouts/{{tissue}}_{{tag}}_layout.txt", zip, tissue=data.tissues, tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/prepMethods/{{tissue}}_{{tag}}_prep_method.txt",zip,tissue=data.tissues,tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam",zip,tissue=data.tissues,tag=data.tags),
        expand(f"{cfg.como_root}/{{tissue}}/geneCounts/{{study}}/{{tissue}}_{{tag}}.tab",zip,tissue=data.tissues,study=data.studies,tag=data.tags),
        expand(f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}/{cfg.sample_filepath.stem}_multiqc_report.html",tissue=set(data.tissues)),
        branch(
            cfg.perform.dump_fastq,
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
                expand(
                    f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_{{end}}_fastqc.zip",
                    zip,
                    tissue=data.tissues_paired,
                    tag=data.tags_paired,
                    end=data.ends_paired
                ),
            ],
            otherwise=[],
        ),
        branch(
            cfg.perform.trim,
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
                expand(
                    f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_{{end}}_fastqc.zip",
                    zip,
                    tissue=data.tissues_paired,
                    tag=data.tags_paired,
                    end=data.ends_paired
                ),
            ],
            otherwise=[],
        ),
        branch(
            cfg.perform.contaminant_screen,
            then=f"{cfg.genome.contaminants_dir}/fastq_screen.conf",
            otherwise=[],
        ),
        branch(
            cfg.perform.fragment_size,
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/fragmentSizes/{{tissue}}_{{tag}}_fragment_size.txt", zip, tissue=data.tissues, tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/fragmentSizes/{{study}}/{{tissue}}_{{tag}}_fragment_size.txt",zip,tissue=data.tissues,study=data.studies,tag=data.tags)
            ],
            otherwise=[]
        ),
        branch(
            cfg.perform.rnaseq_metrics,
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/picard/rnaseq/{{tissue}}_{{tag}}_rnaseq.txt", zip, tissue=data.tissues, tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/strandedness/{{study}}/{{tissue}}_{{tag}}_strandedness.txt",zip,tissue=data.tissues,tag=data.tags,study=data.studies),
            ],
            otherwise=[],
        ),
        branch(
            cfg.perform.insert_size,
            then=[
                expand(f"{cfg.data_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.txt",zip,tissue=data.tissues,tag=data.tags),
                expand(f"{cfg.data_root}/{{tissue}}/picard/hist/{{tissue}}_{{tag}}_insert_size_histo.pdf",zip,tissue=data.tissues,tag=data.tags),
                expand(f"{cfg.como_root}/{{tissue}}/insertSizeMetrics/{{study}}/{{tissue}}_{{tag}}_insert_size.txt",zip,tissue=data.tissues,tag=data.tags,study=data.studies),
            ],
            otherwise=[]
        )

rule copy_config:
    input:
        "config.yaml"
    output:
        f"{cfg.data_root}/{{tissue}}/{{tissue}}_config.yaml"
    resources:
        mem_mb=256,
        runtime=1,
        tissue=""  # intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
    shell: "cp --verbose {input} {output}"

rule preroundup:
    output:
        layout=f"{cfg.data_root}/{{tissue}}/layouts/{{tissue}}_{{tag}}_layout.txt",
        preparation=f"{cfg.data_root}/{{tissue}}/prepMethods/{{tissue}}_{{tag}}_prep_method.txt",
    params:
        sample_name=lambda wildcards: f"{wildcards.tissue}_{wildcards.tag}",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 1 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/preroundup/preroundup_{{tissue}}_{{tag}}.benchmark", cfg.benchmark_count)
    run:
        # example row: SRR12873784,effectorcd8_S1R1,PE,total
        sample_row: pd.Series = data.samples[data.samples["sample"].eq(params.sample_name)]
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
        primary_assembly=f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_{cfg.genome.type}.fa",
        primary_assembly_index=f"{cfg.genome.species_dir}/{cfg.species_name}_{cfg.genome.ensembl_release}_{cfg.genome.type}.fa.fai",
    conda: "envs/generate_genome.yaml"
    threads: 1
    resources:
        mem_mb=8096,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
        network_slots=1
    log: f"{cfg.logs_root}/rule_download_genome_{cfg.species_name}_{cfg.genome.ensembl_release}.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/rule_download_genome_{cfg.species_name}_{cfg.genome.ensembl_release}.benchmark", cfg.benchmark_count)
    shell:
        """
        python3 utils/download_genome.py \
            --taxon-id {cfg.genome.taxon_id} \
            --release-number {cfg.genome.version} \
            --type {cfg.genome.type} \
            --root-save-dir {cfg.genome.species_dir} 1>{log} 2>&1
        """

rule download_contaminant_genomes:
    output:
        done=f"{cfg.genome.contaminants_dir}/.download_complete",
        config=f"{cfg.genome.contaminants_dir}/fastq_screen.conf",
    threads: 1
    params:
        root_output=directory(cfg.genome.contaminants_dir),
    resources:
        mem_mb=6144,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
        network_slots=1
    conda: "envs/screen.yaml"
    log: f"{cfg.logs_root}/download_contaminant_genomes.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/rule_download_contaminant_genomes_{cfg.species_name}.benchmark",cfg.benchmark_count)
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
        mv --verbose "{output.config}.tmp" {output.config} 1>>{log} 2>&1

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
    threads: 10
    resources:
        # +50% of base memory per attempt
        mem_mb=lambda wildcards, attempt: int(51200 * (1 + 0.5 * (attempt - 1))),
        runtime=lambda wildcards, attempt: 150 * attempt,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
    conda: "envs/star.yaml"
    log: f"{cfg.logs_root}/star_index_genome.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/star_index_genome/star_index_genome_{cfg.species_name}.benchmark",cfg.benchmark_count)
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
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tissue="",# intentionally left blank; reference: github.com/jdblischak/smk-simple-slurm/issues/20
    conda: "envs/gffread.yaml"
    log: f"{cfg.logs_root}/generate_transcriptome_fasta.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/rule_generate_transcriptome_{cfg.species_name}.benchmark",cfg.benchmark_count)
    shell: "gffread -w {output.transcriptome} -g {input.genome} {input.gtf} 1>{log} 2>&1"

rule fastq_dump_paired:
    input:
        cfg.sample_filepath
    output:
        r1=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_1.fastq.gz",
        r2=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_2.fastq.gz",
    params:
        srr=lookup(query="sample == '{tissue}_{tag}'", within=data.samples, cols="srr")
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
        network_slots=1
    conda: "envs/SRAtools.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/fastq_dump/{{tissue}}_{{tag}}_fastq_dump.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/fastq_dump_paired/fastq_dump_paired_{{tissue}}_{{tag}}.benchmark", cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        sra_temp="$tmpdir/{wildcards.tissue}_{wildcards.tag}.sra"
        prefetch --max-size u --progress --log-level info --output-file "$sra_temp" {params.srr} 1>{log} 2>&1

        tmp_forward="$tmpdir/{wildcards.tissue}_{wildcards.tag}_1.fastq"
        tmp_reverse="$tmpdir/{wildcards.tissue}_{wildcards.tag}_2.fastq"
        fasterq-dump --force --split-files --progress --threads {threads} --temp "$tmpdir" --outdir "$tmpdir" "$sra_temp" 1>>{log} 2>&1
        pigz --processes {threads} --force "$tmp_forward" "$tmp_reverse"

        printf "\n\n" >> {log}
        mv --verbose "$tmp_forward.gz" "{output.r1}"  1>>{log} 2>&1 &
        mv --verbose "$tmp_reverse.gz" "{output.r2}"  1>>{log} 2>&1 &

        wait
        """

rule fastq_dump_single:
    input:
        cfg.sample_filepath
    output:
        S=f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_S.fastq.gz"
    params:
        srr=lookup(query="sample == '{tissue}_{tag}'",within=data.samples,cols="srr")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    threads: 4
    conda: "envs/SRAtools.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/fastq_dump/{{tissue}}_{{tag}}_fastq_dump.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/fastq_dump_single/fastq_dump_single_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        sra_temp="$tmpdir/{wildcards.tissue}_{wildcards.tag}.sra"
        prefetch --max-size u --progress --log-level info --output-file "$sra_temp" {params.srr} 1>>{log} 2>&1

        tmpfile="$tmpdir/{wildcards.tissue}_{wildcards.tag}.fastq"
        fasterq-dump --force --concatenate-reads --progress --threads {threads} --temp "$tmpdir" --outdir "$tmpdir" "$sra_temp" 1>>{log} 2>&1
        printf "\nGzipping $tmpfile file\n\n" >> {log}
        pigz -6 --processes 4 --force "$tmpfile"

        mv --verbose "$tmpfile.gz" {output.S} 1>>{log} 2>&1
        """


def qc_raw_fastq_paired_input(wildcards):
    if cfg.perform.dump_fastq:
        # We will get data from NCBI, therefore we can get results from fastq_dump_paired directly
        return [rules.fastq_dump_paired.output.r1, rules.fastq_dump_paired.output.r2]

    # Otherwise, we must search the local directory for fiels and return those
    if cfg.local_fastq_filepath:
        sample_name = f"{wildcards.tissue}_{wildcards.tag}"
        for path, subdir, files in os.walk(cfg.local_fastq_filepath):
            for file in files:
                if sample_name in data.sample_names and file.startswith("_".join(wildcards)):
                    forward_read = f"{path}/{file}"
                    reverse_read = forward_read.replace("_1.fastq.gz", "_2.fastq.gz")
                    return [forward_read, reverse_read]
    else:
        raise FileNotFoundError(f"Unable to find directory 'LOCAL_FASTQ_FILES' defined in 'config.yaml'. Attempted searching: {cfg.local_fastq_filepath}")

    print(f"Unable to find files for `qc_raw_fastq_paired` for tissue={wildcards.tissue}, tag={wildcards.tag}")
    return []

rule qc_raw_fastq_paired:
    input:
        reads=qc_raw_fastq_paired_input
    output:
        r1_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_1_fastqc.zip",
        r1_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_1_fastqc.html",
        r2_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_2_fastqc.zip",
        r2_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_2_fastqc.html",
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    threads: 4
    conda: "envs/fastqc.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/raw/{{tissue}}_{{tag}}_fastqc.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/qc_raw_fastq_paired/qc_raw_fastq_paired_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input.reads} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip" "{output.r1_zip}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html" "{output.r1_html}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip" "{output.r2_zip}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html" "{output.r2_html}" 1>>{log} 2>&1
        """

rule qc_raw_fastq_single:
    input:
        rules.fastq_dump_single.output,
    output:
        s_zip=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_S_fastqc.zip",
        s_html=f"{cfg.data_root}/{{tissue}}/fastqc/raw/raw_{{tissue}}_{{tag}}_S_fastqc.html"
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    conda: "envs/fastqc.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/raw/{{tissue}}_{{tag}}_fastqc.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/qc_raw_fastq_single/qc_raw_fastq_single_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads 5 -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.zip" "{output.s_zip}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.html" "{output.s_html}" 1>>{log} 2>&1
        """
def trim_paired_input(wildcards) -> dict[Literal["r1"] | Literal["r2"], str | list[str]]:
    if cfg.perform.dump_fastq:
        return {"r1": rules.fastq_dump_paired.output.r1, "r2": rules.fastq_dump_paired.output.r2}

    if cfg.local_fastq_filepath and cfg.local_fastq_filepath.exists():
        sample_name = f"{wildcards.tissue}_{wildcards.tag}"
        sample_files = sorted(cfg.fastq_files(filter_by=sample_name))
        if len(sample_files) != 2:
            raise ValueError(
                f"Expected 2 FASTQ files for sample '{sample_name}', but found {len(sample_files)} files. "
                f"File(s): {','.join(i.as_posix() for i in sample_files)}"
            )
        return {"r1": sample_files[0].as_posix(),  "r2": sample_files[1].as_posix()}

    print(f"Unable to find any files for `trim_paired` with tissue={wildcards.tissue}, tag={wildcards.tag}")
    return {"r1": [], "r2": []}

rule trim_paired:
    input:
        unpack(trim_paired_input),  # gives 'r1' and 'r2' keywords
    output:
        r1_fastq=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_1.fastq.gz",
        r1_report=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_1_trimming_report.txt",
        r2_fastq=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_2.fastq.gz",
        r2_report=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_2_trimming_report.txt",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 4
    conda: "envs/trim.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_S_trim.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/trim_paired/trim_paired_{{tissue}}_{{tag}}_S.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT
        trim_galore --paired --cores 4 -o "$tmpdir" {input.r1} {input.r2} 1>{log} 2>&1

        printf "\n\n" >> {log}
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_val_1.fq.gz" "{output.r1_fastq}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1.fastq.gz_trimming_report.txt" "{output.r1_report}" 1>>{log} 2>&1

        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_val_2.fq.gz" "{output.r2_fastq}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2.fastq.gz_trimming_report.txt" "{output.r2_report}" 1>>{log} 2>&1
        """


def trim_single_input(wildcards) -> list[str]:
    if cfg.perform.dump_fastq:
        return [rules.fastq_dump_single.output.S]

    if cfg.local_fastq_filepath and cfg.local_fastq_filepath.exists():
        sample_name = f"{wildcards.tissue}_{wildcards.tag}"
        sample_files = cfg.fastq_files(filter_by=sample_name)
        if len(sample_files) != 1:
            raise ValueError(
                f"Expected 1 FASTQ files for sample '{sample_name}', but found {len(sample_files)} files. "
                f"File(s): {','.join(i.as_posix() for i in sample_files)}"
            )
        return [sample_files[0].as_posix()]

    print(f"Unable to find any files for `trim_single` with tissue={wildcards.tissue}, tag={wildcards.tag}")
    return []

rule trim_single:
    input:
        S=trim_single_input
    output:
        S_fastq=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_S.fastq.gz",
        S_report=f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_S_trimming_report.txt"
    # See the trim_galore `--cores` setting for details on why 16 was chosen
    # https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 4
    conda: "envs/trim.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_trim.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/trim_single/trim_single_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir="$(mktemp -d)"
        trap "rm -rf $tmpdir" EXIT

        trim_galore --cores 4 --output_dir "$tmpdir" {input.S} 1>{log} 2>&1

        mkdir --verbose -p "$(dirname {output.S_fastq})" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_trimmed.fq.gz" "{output.S_fastq}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_S.fastq.gz_trimming_report.txt" "{output.S_report}" 1>>{log} 2>&1
        """

rule qc_trim_fastq_paired:
    input:
        r1=rules.trim_paired.output.r1_fastq,
        r2=rules.trim_paired.output.r2_fastq
    output:
        r1_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_1_fastqc.zip",
        r1_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_1_fastqc.html",
        r2_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_2_fastqc.zip",
        r2_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_2_fastqc.html",
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    conda: "envs/trim.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/trimmed/{{tissue}}_{{tag}}_fastqc.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/qc_trim_fastq_paired/qc_trim_fastq_paired_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.zip" "{output.r1_zip}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_1_fastqc.html" "{output.r1_html}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.zip" "{output.r2_zip}" 1>>{log} 2>&1
        mv --verbose "$tmpdir/{wildcards.tissue}_{wildcards.tag}_2_fastqc.html" "{output.r2_html}" 1>>{log} 2>&1
        """

rule qc_trim_fastq_single:
    input:
        rules.trim_single.output.S_fastq,
    output:
        s_zip=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_S_fastqc.zip",
        s_html=f"{cfg.data_root}/{{tissue}}/fastqc/trimmed/trimmed_{{tissue}}_{{tag}}_S_fastqc.html"
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    conda: "envs/trim.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fastqc/trimmed/{{tissue}}_{{tag}}_fastqc.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/qc_trim_fastq_single/qc_trim_fastq_single_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        r"""
        tmpdir=$(mktemp -d)
        tmp_zip="$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.zip"
        tmp_html="$tmpdir/{wildcards.tissue}_{wildcards.tag}_S_fastqc.html"
        trap "rm -rf $tmpdir" EXIT

        fastqc {input} --threads {threads} -o "$tmpdir" 1>{log} 2>&1

        printf "\n\n" >> {log}
        mv --verbose "$tmp_zip" "{output.s_zip}" 1>>{log} 2>&1
        mv --verbose "$tmp_html" "{output.s_html}" 1>>{log} 2>&1
        """

def align_input(wildcards):
    sample_name = f"{wildcards.tissue}_{wildcards.tag}"

    # find local files
    if not cfg.perform.trim and not cfg.perform.dump_fastq:
        return [i.as_posix() for i in cfg.local_fastq_filepath.rglob(f"{sample_name}_[12S].fastq.gz")]

    _endtype: str = data.samples.loc[data.samples["sample"] == sample_name, "endtype"].values[0]
    # Get files from trim_paired and trim_single
    if cfg.perform.trim:
        if _endtype == "PE":
            return expand(rules.trim_paired.output.r1_fastq, **wildcards) + expand(rules.trim_paired.output.r2_fastq, **wildcards)
        elif _endtype == "SE":
            return expand(rules.trim_single.output.S_fastq, **wildcards)
        else:
            raise ValueError(f"Invalid endtype '{_endtype}' for sample '{sample_name}'. Must be one of 'PE' or 'SE'.")

    # get files from dump_fastq_paired and dump_fastq_single
    if cfg.perform.dump_fastq:
        if _endtype == "PE":
            return expand(rules.fastq_dump_paired.output.r1, **wildcards) + expand(rules.fastq_dump_paired.output.r2, **wildcards)
        elif _endtype == "SE":
            return expand(rules.fastq_dump_single.output.S, **wildcards)
        else:
            raise ValueError(f"Invalid endtype '{_endtype}' for sample '{sample_name}'. Must be one of 'PE' or 'SE'.")

    print(f"Unable to find any files for `align` with tissue={wildcards.tissue}, tag={wildcards.tag}")
    return []


rule align:
    input:
        files=align_input,
        genome=rules.download_genome.output.primary_assembly
    output:  # Not all outputs listed are used, but they are listed so Snakemake knows to clean them up on workflow reruns
        bam_file=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam",
        final_log=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.final.out",
        intermediate_log=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.out",
        progress=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Log.progress.out",
        splice_junctions=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_SJ.out.tab",
        gene_table=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.tab",
        tx_bam=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.toTranscriptome.bam"
    params:
        star_genome=f"{cfg.genome.species_dir.as_posix().removesuffix('/')}/star",
        gene_table=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_ReadsPerGene.out.tab",
        bam_output=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Aligned.sortedByCoord.out.bam",
        tx_bam_output=f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}_Aligned.toTranscriptome.out.bam",
    resources:
        mem_mb=lambda wildcards, attempt: 40960 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 5
    conda: "envs/star.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/align/{{tissue}}_{{tag}}_star_align.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/align/align_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
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
        mv --verbose $tmpdir/* "$(dirname {output.gene_table})/" 1>>{log} 2>&1
        mv --verbose {params.gene_table} {output.gene_table} 1>>{log} 2>&1
        mv --verbose {params.bam_output} {output.bam_file} 1>>{log} 2>&1
        mv --verbose {params.tx_bam_output} {output.tx_bam} 1>>{log} 2>&1
        """

rule index_bam_file:
    input:
        rules.align.output.bam_file,
    output:
        f"{cfg.data_root}/{{tissue}}/align/{{tag}}/{{tissue}}_{{tag}}.bam.bai",
    resources:
        mem_mb=lambda wildcards, attempt: 512 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 4
    conda: "envs/samtools.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/align/{{tissue}}_{{tag}}_index_bam.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/index_bam_file/index_bam_filepath_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
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
    resources:
        mem_mb=lambda wildcards, attempt: 32768 * attempt,
        runtime=lambda wildcards, attempt: 40 * attempt,
        tissue=lambda wildcards: wildcards.tissue
    threads: 8
    conda: "envs/salmon.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/salmon_quant/{{tissue}}_{{tag}}_salmon_quant.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/salmon_quantification/salmon_quantification_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
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
        mv --verbose {params.outdir}/quant.sf {output.quant} 1>>{log} 2>&1
        mv --verbose {params.outdir}/cmd_info.json {output.meta} 1>>{log} 2>&1
        """


def contaminant_screen_input_paired(wildcards) -> dict[Literal["screen_config"] | Literal["files"], list[str]]:
    returns: dict[Literal["screen_config"] | Literal["files"], list[str]] = {"screen_config": [], "files": []}
    if not cfg.perform.contaminant_screen:
        return returns

    sample_name = f"{wildcards.tissue}_{wildcards.tag}"
    pe_samples: pd.DataFrame = data.samples[(data.samples["sample"] == sample_name) & (data.samples["endtype"] == "PE")]
    if pe_samples.empty:
        return returns

    if cfg.perform.trim:
        returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
        returns["files"] = [rules.trim_paired.output.r1_fastq, rules.trim_paired.output.r2_fastq]
    elif cfg.perform.dump_fastq:
        returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
        returns["files"] = [rules.fastq_dump_paired.output.r1, rules.fastq_dump_paired.output.r2]
    else:
        files = [i.as_posix() for i in cfg.local_fastq_filepath.rglob(f"{sample_name}_[12].fastq.gz")]
        if files:
            returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
            returns["files"] = files
        else:
            raise FileNotFoundError(f"Unable to find paired-end FASTQ files for sample '{sample_name}' in directory '{cfg.local_fastq_filepath}'")

    return returns

rule contaminant_screen_paired:
    input:
        unpack(contaminant_screen_input_paired)
    output:
        r1=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_1_screen.txt",
        r2=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_2_screen.txt"
    params:
        output_dir=f"{cfg.data_root}/{{tissue}}/fq_screen",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 5
    conda: "envs/screen.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_fastq_screen.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/contaminant_screen_paired/contaminant_screen_paired_{{tissue}}_{{tag}}.benchmark", cfg.benchmark_count)
    shell:
        """
        outdir=$(dirname {output.r1})
        mkdir -p "$outdir"
        fastq_screen --force --aligner Bowtie2 --threads {threads} --conf {input.screen_config} --outdir "$outdir" {input.files} 1>{log} 2>&1
        """


def contaminant_screen_input_single(wildcards) -> dict[Literal["screen_config"] | Literal["files"], list[str]]:
    returns: dict[Literal["screen_config"] | Literal["files"], list[str]] = {"screen_config": [], "files": []}
    if not cfg.perform.contaminant_screen:
        return returns

    sample_name = f"{wildcards.tissue}_{wildcards.tag}"
    se_samples: pd.DataFrame = data.samples[(data.samples["sample"] == sample_name) & (data.samples["endtype"] == "SE")]
    if se_samples.empty:
        return returns

    if cfg.perform.trim:
        returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
        returns["files"] = [rules.trim_single.output.S_fastq]
        return returns
    elif cfg.perform.dump_fastq:
        returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
        returns["files"] = [rules.fastq_dump_single.output.S]
        return returns
    else:
        files = [i.as_posix() for i in cfg.local_fastq_filepath.rglob(f"{sample_name}_S.fastq.gz")]
        if files:
            returns["screen_config"] = [rules.download_contaminant_genomes.output.config]
            returns["files"] = files
            return returns
        else:
            raise FileNotFoundError(f"Unable to find single-end FASTQ file for sample '{sample_name}' in directory '{cfg.local_fastq_filepath}'")

    print(f"Unable to find any files for `contaminant_screen_single` with tissue={wildcards.tissue}, tag={wildcards.tag}")
    return returns

rule contaminant_screen_single:
    input:
        unpack(contaminant_screen_input_single)  # gives `input.screen_config` and `input.files`
    output:
        S=f"{cfg.data_root}/{{tissue}}/fq_screen/{{tissue}}_{{tag}}_S_screen.txt"
    params:
        output_dir=f"{cfg.data_root}/{{tissue}}/fq_screen"
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 20 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 5
    conda: "envs/screen.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/contaminant_screen/{{tissue}}_{{tag}}_fastq_screen.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/contaminant_screen_single/contaminant_screen_single{{tissue}}_{{tag}}.benchmark", cfg.benchmark_count)
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
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        runtime=lambda wildcards, attempt: 20 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    conda: "envs/rseqc.yaml"
    threads: 4
    log: f"{cfg.logs_root}/{{tissue}}/fragment_size/{{tissue}}_{{tag}}_fragment_size.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/fragment_size/fragment_size_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell: "python3 utils/get_fragment_size.py --input {input.bam} --bai {input.bai} --refgene {input.bed_file} --log {log} --output {output}"


rule insert_size:
    input:
        bam=rules.align.output.bam_file,
        layout=rules.preroundup.output.layout
    output:
        txt=f"{cfg.data_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.txt",
        pdf=f"{cfg.data_root}/{{tissue}}/picard/hist/{{tissue}}_{{tag}}_insert_size_histo.pdf",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 1
    conda: "envs/picard.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/picard/insert/{{tissue}}_{{tag}}_insert_size.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/insert_size/insert_size_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
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
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    conda: "envs/picard.yaml"
    threads: 1
    log: f"{cfg.logs_root}/{{tissue}}/picard/rnaseq/{{tissue}}_{{tag}}_rnaseq_metrics.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/rnaseq_metrics/rnaseq_metrics_{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
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
    input:
        rules.fragment_size.output,
    output:
        f"{cfg.como_root}/{{tissue}}/fragmentSizes/{{sample}}/{{tissue}}_{{tag}}_fragment_size.txt"
    resources:
        mem_mb=256,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/copy_fragment_size/copy_fragment_size_{{sample}}/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """cp --verbose {input} {output}"""

rule copy_insert_size:
    input:
        rules.insert_size.output.txt,
    output:
        f"{cfg.como_root}/{{tissue}}/insertSizeMetrics/{{sample}}/{{tissue}}_{{tag}}_insert_size.txt"
    resources:
        mem_mb=256,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/copy_insert_size/copy_insert_size_{{sample}}/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """cp --verbose {input} {output}"""

rule copy_rnaseq_metrics:
    input:
        rules.rnaseq_metrics.output.strand,
    output:
        f"{cfg.como_root}/{{tissue}}/strandedness/{{sample}}/{{tissue}}_{{tag}}_strandedness.txt"
    resources:
        mem_mb=256,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/copy_rnaseq_metrics/copy_rnaseq_metrics_{{sample}}/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """cp --verbose {input} {output}"""


rule copy_gene_counts:
    input:
        rules.align.output.gene_table,
    output:
        f"{cfg.como_root}/{{tissue}}/geneCounts/{{sample}}/{{tissue}}_{{tag}}.tab"
    resources:
        mem_mb=256,
        runtime=1,
        tissue=lambda wildcards: wildcards.tissue,
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/copy_gene_counts/copy_gene_counts_{{sample}}/{{tissue}}_{{tag}}.benchmark",cfg.benchmark_count)
    shell:
        """cp --verbose {input} {output}"""


def multiqc_contamination_input(wildcards) -> list[str]:
    if not cfg.perform.contaminant_screen:
        return []

    pe_samples = data.samples.loc[(data.samples["sample"].str.contains(wildcards.tissue)) & (data.samples["endtype"] == "PE")]
    se_samples = data.samples.loc[(data.samples["sample"].str.contains(wildcards.tissue)) & (data.samples["endtype"] == "SE")]
    files: list[str] = []
    if not pe_samples.empty:
        tissues, tags = pe_samples["sample"].str.split(n=1, pat="_", expand=True).T.values
        files += expand(rules.contaminant_screen_paired.output.r1, zip, tissue=tissues, tag=tags)
        files += expand(rules.contaminant_screen_paired.output.r2, zip, tissue=tissues, tag=tags)
    if not se_samples.empty:
        tissues, tags = se_samples["sample"].str.split(n=1, pat="_", expand=True).T.values
        files += expand(rules.contaminant_screen_single.output.S, zip, tissue=tissues, tag=tags)
    return files


rule multiqc:
    input:
        raw_fastq=lambda wildcards: [] if not cfg.perform.dump_fastq else expand(f"{cfg.data_root}/{{tissue}}/raw/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
        trimmed_fastq=lambda wildcards: [] if not cfg.perform.trim else expand(f"{cfg.data_root}/{{tissue}}/trim/{{tissue}}_{{tag}}_{{end}}.fastq.gz",zip,tissue=data.tissues_paired,tag=data.tags_paired,end=data.ends_paired),
        aligned_fastq=expand(rules.align.output.bam_file, zip, tissue=data.tissues, tag=data.tags),
        contaminantion=multiqc_contamination_input,
        insert_sizes=lambda wildcards: [] if not cfg.perform.insert_size else expand(rules.insert_size.output.txt,zip,tissue=data.tissues,tag=data.tags),
        rnaseq_metrics=lambda wildcards: [] if not cfg.perform.rnaseq_metrics else expand(rules.rnaseq_metrics.output.metrics, zip, tissue=data.tissues, tag=data.tags),
        fragment_sizes=lambda wildcards: [] if not cfg.perform.fragment_size else expand(rules.fragment_size.output, zip, tissue=data.tissues, tag=data.tags),
        salmon_quant=expand(rules.salmon_quantification.output.quant, zip, tissue=data.tissues, tag=data.tags),
    output:
        output_file=f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}/{cfg.sample_filepath.stem}_multiqc_report.html",
    params:
        config_file_basename=cfg.sample_filepath.stem,
        input_directory=f"{cfg.data_root}/{{tissue}}",
        output_directory=directory(f"{cfg.data_root}/{{tissue}}/multiqc/{cfg.sample_filepath.stem}"),
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tissue=lambda wildcards: wildcards.tissue,
    threads: 1
    conda: "envs/multiqc.yaml"
    log: f"{cfg.logs_root}/{{tissue}}/multiqc/{{tissue}}_multiqc.log"
    benchmark: repeat(f"{cfg.benchmark_dir}/{{tissue}}/multiqc/multiqc_{{tissue}}.benchmark",cfg.benchmark_count)
    shell:
        """
        filename="$(basename {output})"
        mkdir -p "{params.output_directory}"

        multiqc --interactive --force --title "{wildcards.tissue}" \
          --filename {params.config_file_basename}_multiqc_report.html \
          --outdir "{params.output_directory}" {params.input_directory} 1>{log} 2>&1
        """
