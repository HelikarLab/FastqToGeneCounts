# ruff: noqa: S310 S321 T201
import argparse
import csv
import ftplib
import gzip
import http
import io
import json
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Literal

_ncbi_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
_ensembl_url = "ftp.ensembl.org"

_ucsc_url = "https://api.genome.ucsc.edu"
_ref_flat_url = "https://hgdownload.soe.ucsc.edu/goldenPath/{final_genome}/database/refFlat.txt.gz"  # fmt: skip

MUS_MUSCULUS: int = 10090
HOMO_SAPIENS: int = 9606
MACACA_MULATTA: int = 9544
MIN_ENSEMBL_RELEASE: int = 19


def _get_cache_dir() -> Path:
    curr_dir = Path(__file__)
    parent = curr_dir.parent
    while "Snakefile" not in [i.name for i in parent.iterdir()]:
        parent = parent.parent
        if parent.as_posix() == "/":
            raise FileNotFoundError(f"Unable to find Snakemake file. Search started in: {curr_dir}")
    return parent


def _get_species_taxon_cache_file() -> Path:
    cache_file = _get_cache_dir() / ".taxon_to_species.json"
    cache_file.touch(exist_ok=True)
    return cache_file


def get_latest_release(*, print_warning: bool = True) -> int:
    pub_root: ftplib.FTP = ftplib.FTP(_ensembl_url)
    pub_root.login()
    latest_release: int = -1
    release: str
    for release in pub_root.nlst("/pub"):
        if release.startswith("/pub/release-") and not release.endswith(".txt"):
            release_number: int = int(release.split("-")[-1])
            latest_release = max(latest_release, release_number)

    retries = 3
    while not pub_root.nlst(f"/pub/release-{latest_release}/README") and retries > 0:
        if print_warning:
            print(f"❗❗❗ WARNING ❗❗❗ Could not find files for release 'release-{latest_release}', moving to 'release-{latest_release - 1}'")
        latest_release -= 1
        retries -= 1
    return latest_release


def species_from_taxon(
    taxon_id: int,
    replace_with: str = "_",
    *,
    lowercase: bool = True,
    replace_spaces: bool = True,
) -> str:
    """Get the species name from the provided taxon ID.

    It assumes the taxon ID is valid, as defined by `_is_valid_taxon`
    It will first read from a JSON cache that is located in the same directory as the Snakemake file.

    Args:
        taxon_id: The taxon ID to convert to species name
        lowercase: If True, the species name will be returned in lowercase
        replace_spaces: If True, spaces will be replaced with underscores
        replace_with: The character to replace spaces with

    Returns:
        str: The species name that corresponds to the ensembl species. For example: 9606 -> homo_sapiens
    """
    cache_file = _get_species_taxon_cache_file()
    cache_dict = {}
    if cache_file.exists():
        text = cache_file.read_text().rstrip("\n")
        cache_dict = json.loads(text) if text else {}

    key = str(taxon_id)
    if key in cache_dict:
        species = cache_dict[key].replace("_", " ").title()
        species = species.replace(" ", replace_with) if replace_spaces else species
        return species.lower() if lowercase else species

    invalid_taxon_error = f"The taxonomy name you specified ({taxon_id}) is not a recognized NCBI Taxonomy name."
    backoff = 2.0
    for attempt in range(6):
        try:
            req = urllib.request.Request(f"{_ncbi_url}/taxonomy/taxon/{taxon_id}/name_report")  # noqa: S310
            with urllib.request.urlopen(req, timeout=30) as response:  # noqa: S310
                as_json = json.loads(response.read().decode("utf-8"))

            # If we get an invalid taxon, set the species to a null/empty value
            is_valid = True
            for response_list in as_json["reports"]:
                if "errors" in response_list:
                    error_keys, error_values = response_list["errors"][0].items()
                    if ("reason" in error_keys) and (invalid_taxon_error in error_values):
                        is_valid = False

            species = ""
            if is_valid:
                species = as_json["reports"][0]["taxonomy"]["current_scientific_name"]["name"]
                species = species.lower() if lowercase else species
                species = species.replace(" ", replace_with) if replace_spaces else species

            cache_dict[key] = species
            cache_file.write_text(json.dumps(cache_dict, indent=2, sort_keys=True))
        except urllib.error.HTTPError as e:  # noqa: PERF203
            status = http.HTTPStatus(e.code)
            if status.is_client_error and attempt < 5:  # noqa: PLR2004
                time.sleep(backoff)
                backoff *= 2
                continue
            raise
        except urllib.error.URLError:
            raise
        else:
            return species
    raise RuntimeError("Failed to resolve species after 5 retries.")


def is_valid_release_number(release_number: str | int) -> bool:
    cache_file = _get_cache_dir() / ".valid_release_numbers.txt"
    cache_file.touch(exist_ok=True)
    releases: list[str] = cache_file.read_text().rstrip("\n").split("\n")
    if "" in releases:
        releases.remove("")
    try:
        if isinstance(release_number, str):
            release_number = int(release_number.split("-")[-1]) if release_number.startswith("release-") else int(release_number)
    except ValueError:  # Unable to convert to an integer
        return False

    if str(release_number) in releases:
        return True

    latest_release = str(get_latest_release(print_warning=False))
    if latest_release not in releases:
        releases.append(latest_release)
    sorted_releaes = sorted(int(i) for i in releases)
    cache_file.write_text("\n".join(str(i) for i in sorted_releaes))
    return True


class NCBI:
    def __init__(
        self,
        taxon_id: int,
        release_number: str,
        fasta_type: Literal["primary_assembly", "toplevel"],
        *,
        show_progress: bool,
    ) -> None:
        """Initialize an NCBI connection.

        :param taxon_id: The taxon ID to download the genome for
        :param release_number: The Ensembl release number to download the genome for
        :param show_progress: Whether to show a progress bar during downloads
        """
        self._taxon_id: int = taxon_id
        self._release_number: str = release_number
        self._species_name: str = species_from_taxon(self._taxon_id)
        self._progress: bool = show_progress
        self._type: Literal["primary_assembly", "toplevel"] = fasta_type

        if not self._release_number.startswith(("latest", "release-")) and not self._release_number.isdigit():
            raise ValueError("""release_number must be an integer or start with 'latest' or 'release-'.\nExamples: "115", "latest", "release-112""")

        taxon_valid = species_from_taxon(taxon_id=taxon_id)  # we will get an empty string or the species name
        if not taxon_valid:
            raise ValueError(
                f"The provided taxon_id ({self._taxon_id}) is not valid.\nPlease validate your taxon ID at https://www.ncbi.nlm.nih.gov/Taxonomy"
            )

        if self._release_number != "latest" and not is_valid_release_number(self._release_number):
            raise ValueError(
                f"The provided release number ({self._release_number}) is not valid. It should be in the format of `release-112` or `112`.\n"
                f"Please validate your release number at https://ftp.ensembl.org/pub"
            )

        # Inputs are valid, we can now download the fasta primary assembly
        if self._release_number == "latest":
            self._release_number = f"release-{get_latest_release()}"
        elif self._release_number.isdigit():
            self._release_number = f"release-{self._release_number}"

        self._ftp = ftplib.FTP(_ensembl_url)
        self._ftp.login()

    def __enter__(self):
        """Enable use of `with` statement."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Enable use of `with` statement."""
        self._ftp.quit()

    def download_primary_assembly(self, save_directory: str) -> None:
        """Download the reference genome from the UCSC Genome Browser and save it as a FASTA file.

        # The full input path of the genome fasta file
        GENOME_FASTA_FILE: "genome/Mus_musculus.GRCm39.dna.primary_assembly.fa"
        GTF_FILE: "genome/Mus_musculus.GRCm39.111.gtf"
        """
        save_directory: Path = Path(save_directory)
        final_output_filepath = Path(save_directory, f"{self._species_name}_{self._release_number}_{self._type}.fa")
        if final_output_filepath.exists():
            return

        save_directory.mkdir(exist_ok=True, parents=True)

        print(f"Downloading FASTA file to: {final_output_filepath}")
        fasta_root = f"/pub/{self._release_number}/fasta/{self._species_name}/dna/"
        assembly_suffix = f".dna.{self._type}.fa.gz"
        assembly_path = ""
        total_download_size: int = 0

        print(f"Searching: {fasta_root}")
        print(f"{assembly_suffix=}")
        for filename in self._ftp.nlst(fasta_root):
            print(f"\tFile: {filename}")
            if filename.endswith(assembly_suffix):
                assembly_path = filename

                size = self._ftp.size(filename)
                if size is not None:
                    total_download_size += size
                break

        # If the primary assembly path is not empty, we can download it
        if assembly_path != "":
            assembly_filename = assembly_path.split("/")[-1]
            with Path(save_directory, assembly_filename).open("wb") as fasta_file:
                self._ftp.retrbinary(f"RETR {assembly_path}", fasta_file.write)
        # Otherwise we must download all *.dna.chromosome.*.fa.gz files and concatenate them
        else:
            assembly_filename = ""
            chromosome_files = []
            for filename in self._ftp.nlst(fasta_root):
                if f".dna.{self._type}.fa.gz" in filename:
                    chromosome_files.append(filename)
                if filename.endswith(f".dna.{self._type}.1.fa.gz"):
                    assembly_filename = filename.split("/")[-1].replace(".1.fa.gz", ".fa.gz")

            output_fasta_file = Path(save_directory, assembly_filename)
            print(f"Output fasta file will be: {output_fasta_file}")
            with output_fasta_file.open("wb") as fasta_file:
                for remote_chromosome in chromosome_files:
                    chromosome_filename = remote_chromosome.split("/")[-1]
                    local_chromosome = Path(save_directory, chromosome_filename)
                    with local_chromosome.open("wb") as chr_out:
                        self._ftp.retrbinary(f"RETR {remote_chromosome}", chr_out.write)

                    with local_chromosome.open("rb") as chr_in:
                        fasta_file.write(chr_in.read())

            # Remove primary assembly files that are not the "primary_assembly"
            for filename in save_directory.iterdir():
                # Skip files not related to primary_assembly
                if self._type not in filename.name:
                    continue
                if not filename.name.endswith(f".dna.{self._type}.fa.gz"):
                    Path(save_directory, filename).unlink(missing_ok=True)

        # Once the primary assembly is downloaded/created, un-gzip it
        txt_filename = assembly_filename.replace(".gz", "")
        with gzip.open(f"{save_directory}/{assembly_filename}", "rb") as f_in, Path(save_directory, txt_filename).open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Remove the gzipped file then convert file to lowercase name
        Path(save_directory, assembly_filename).unlink()
        shutil.move(Path(save_directory, txt_filename), final_output_filepath)

    def download_gtf_file(self, save_directory: str) -> None:
        """Download a GTF file.

        :param: save_directory: The directory to save the GTF file to.

        """
        save_directory: Path = Path(save_directory)
        final_output_filepath = Path(save_directory, f"{self._species_name}_{self._release_number}.gtf")
        if final_output_filepath.exists():
            return
        print(f"Downloading GTF file to: {final_output_filepath}")
        save_directory.mkdir(parents=True, exist_ok=True)

        gtf_root = f"/pub/{self._release_number}/gtf/{self._species_name}"
        gtf_suffix = f"{self._release_number.split('-')[-1]}.gtf.gz"

        # Download file from ensembl FTP server
        # TODO: add a progress bar to this download
        gtf_filename = ""
        for filename in self._ftp.nlst(gtf_root):
            if filename.endswith(gtf_suffix):
                gtf_filename = filename.split("/")[-1]
                with Path(save_directory, gtf_filename).open("wb") as gtf_file:
                    self._ftp.retrbinary(f"RETR {filename}", gtf_file.write)
                break

        # ungzip the gtf file
        txt_filename = gtf_filename.replace(".gz", "")
        txt_obj = Path(save_directory, txt_filename)
        with gzip.open(f"{save_directory}/{gtf_filename}", "rb") as f_in:
            with txt_obj.open("wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove the gzipped file then convert file to lowercase name
        Path(save_directory, gtf_filename).unlink()
        shutil.move(txt_obj, final_output_filepath)

    def download_cdna_fasta_file(self, save_directory: str) -> None:
        """Download Ensembl cDNA FASTA (transcriptome) for a species + release."""
        save_directory: Path = Path(save_directory)
        save_directory.mkdir(parents=True, exist_ok=True)

        final_output = Path(save_directory, f"{self._species_name}_{self._release_number}_cdna.fa")
        if final_output.exists():
            return

        cdna_root = f"/pub/{self._release_number}/fasta/{self._species_name}/cdna"

        # Prefer *.cdna.all.fa.gz (standard transcriptome)
        cdna_gz_remote = None
        for remote in self._ftp.nlst(cdna_root):
            if remote.endswith(".cdna.all.fa.gz"):
                cdna_gz_remote = remote
                break
        if cdna_gz_remote is None:
            raise FileNotFoundError(
                f"Could not find a .cdna.all.fa.gz file under {cdna_root} for {self._species_name} and release {self._release_number}"
            )

        gz_name = cdna_gz_remote.split("/")[-1]
        local_gz = Path(save_directory, gz_name)
        with local_gz.open("wb") as out:
            self._ftp.retrbinary(f"RETR {cdna_gz_remote}", out.write)

        # Un-zip to final_output
        with gzip.open(local_gz, "rb") as f_in, final_output.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        local_gz.unlink(missing_ok=True)


def ref_flat_file_creation(taxon_id: int, save_directory: str) -> None:
    """Download the `refFlat.txt` file from the UCSC Genome Browser and save it as a file.

    :param taxon_id: The taxon ID to download the `refFlat.txt` file for
    :param save_directory: The directory to save the `refFlat.txt` file
    """
    species_name_with_space = species_from_taxon(taxon_id, replace_spaces=False)
    species_name = species_name_with_space.replace(" ", "_")
    save_directory: Path = Path(save_directory)
    save_directory.mkdir(parents=True, exist_ok=True)

    final_output_filepath = Path(save_directory, f"{species_name}_ref_flat.txt")
    if final_output_filepath.exists():
        return

    print(f"Saving Reference FLAT file to: {final_output_filepath}")
    # Get a list of genomes
    req = urllib.request.Request(f"{_ucsc_url}/list/ucscGenomes")
    with urllib.request.urlopen(req, timeout=30) as response:
        as_json = json.loads(response.read())

    genome_prefix: str = ""
    genome_version: int = -1

    for genome_id in as_json["ucscGenomes"]:
        genome_data: dict = as_json["ucscGenomes"][genome_id]
        ucsc_species_name = genome_data["scientificName"]

        if ucsc_species_name.lower() == species_name_with_space:
            prefix = re.search(r"^[a-zA-Z]+", genome_id)
            version = re.search(r"\d+", genome_id)

            if version is not None and prefix is not None:
                version = int(version.group())
                if version > genome_version:
                    genome_version = version
                    genome_prefix = prefix.group()

    if genome_version == -1 and genome_prefix == "":
        raise ValueError(f"Could not find UCSC genome for taxon ID: {taxon_id}")

    # Download `ref_flat_url` and save it as a file without creating an intermeidate gzip file
    final_genome: str = f"{genome_prefix}{genome_version}"
    download_url = _ref_flat_url.format(final_genome=final_genome)
    txt_file = Path(save_directory, f"ref_flat_{final_genome}.txt")

    # TODO: add a progress bar for this download
    try:
        response = urllib.request.urlopen(download_url)
    except urllib.error.HTTPError:
        print(f"Unable to find url: '{download_url}'")
        raise

    with gzip.GzipFile(fileobj=io.BytesIO(response.read()), mode="rb") as decompressed, txt_file.open("wb") as o_stream:
        o_stream.write(decompressed.read())

    shutil.move(txt_file, final_output_filepath)


def r_rna_interval_list_creation(taxon_id: int, fasta_type: Literal["primary_assembly", "toplevel"], save_directory: str) -> None:
    """Create an interval file list suitable for CollectRnaSeqMetrics.

    It should be executed at the end of the genome generation process so the required `.fa` file is availale

    Initially created by Kamil Slowikowski (December 12, 2014)
    Modified by: Arindam Ghosh (July 24, 2019)
    Modified by: Josh Loecker (April 13, 2023)

    :param taxon_id: The taxon ID to create the rRNA interval list for
    :param fasta_type: The type of FASTA file to use (either "primary_assembly" or "toplevel")
    :param save_directory: The directory to save the rRNA interval list
    :raises FileNotFoundError: «The primary assembly file could not be found» if the primary assembly file is not found
    """
    species_name = species_from_taxon(taxon_id)
    save_directory: Path = Path(save_directory)
    save_directory.mkdir(parents=True, exist_ok=True)

    final_output_filepath = Path(save_directory, f"{species_name}_rrna.interval_list")
    if final_output_filepath.exists():
        return

    print(f"Saving rRNA interval list to: {final_output_filepath}")

    genome_sizes = Path(save_directory, f"{species_name}_genome_sizes.txt")
    genes: str = next(file.name for file in save_directory.iterdir() if file.name.endswith(".gtf"))
    primary_assembly_fa: str = next(file.name for file in save_directory.iterdir() if file.name.endswith(f"{fasta_type}.fa"))

    genes: Path = Path(save_directory, genes)
    primary_assembly_fa: Path = Path(save_directory, primary_assembly_fa)
    primary_assembly_fai: Path = Path(f"{primary_assembly_fa}.fai")

    # If primary_assembly_fa doesn't exist, show error and quit
    if not primary_assembly_fa.exists():
        raise FileNotFoundError(f"The primary assembly file could not be found\nSearching for: '{primary_assembly_fa}'")  # fmt: skip

    subprocess.run(  # noqa: S603
        ["/usr/bin/env", "samtools", "faidx", primary_assembly_fa],
        capture_output=True,
        check=False,
        shell=False,
    )

    with primary_assembly_fai.open(mode="r") as primary_assembly_fai_in, genome_sizes.open(mode="w") as genome_sizes_out:
        for line in primary_assembly_fai_in:
            fields = line.split()
            genome_sizes_out.write(f"{fields[0]}\t{fields[1]}\n")

    r_rna_interval_list = Path(save_directory, f"{species_name}_rrna.interval_list")
    with r_rna_interval_list.open(mode="w") as o_stream:
        with genome_sizes.open(mode="r") as i_stream:
            for line in i_stream:
                fields = line.split()
                o_stream.write(f"@SQ\tSN:{fields[0]}\tLN:{fields[1]}\tUR:file:{primary_assembly_fa}\n")  # fmt: skip

        with genes.open(mode="r") as i_stream:
            for line in i_stream:
                if 'gene_biotype "rRNA"' in line:
                    fields = line.split()
                    if fields[2] == "gene":
                        gene_id = fields[-1].split('"')[1]
                        o_stream.write(f"{fields[0]}\t{fields[3]}\t{fields[4]}\t{fields[6]}\t{gene_id}\n")

    shutil.move(r_rna_interval_list, final_output_filepath)

    # Sort the interval list
    # TODO: Is this needed? As of now, I don't think so
    # subprocess.run(["sort", "-k1V", "-k2n", "-k3n", "-o", r_rna_interval_list, r_rna_interval_list])


def bed_file_creation(taxon_id: int, save_directory: str) -> None:
    """Convert USCS refFlat.txt files into BED format.

    # The reference BED file built from the GTF for RSEQC option
    BED_FILE: "genome/mm10_RefSeq.bed"

    modified from: informationsea at https://gist.github.com/informationsea/439d4fc53ea2b17cfb05

    :param taxon_id: The taxon ID to create the BED file for
    :param save_directory: The directory to save the BED file
    """
    species_name = species_from_taxon(taxon_id)
    save_directory: Path = Path(save_directory)
    bed_output_file = Path(save_directory, f"{species_name}.bed")
    save_directory.mkdir(parents=True, exist_ok=True)

    if bed_output_file.exists():
        return

    print(f"Saving BED file to: {bed_output_file}")

    ref_flat_file = Path(save_directory, f"{species_name}_ref_flat.txt")
    if not ref_flat_file.exists():
        raise FileNotFoundError(f"The refFlat file could not be found\nSearching for: '{ref_flat_file}'")

    reader = csv.reader(ref_flat_file.open(mode="r"), delimiter="\t", quotechar=None)
    writer = csv.writer(bed_output_file.open(mode="w"), delimiter="\t", quotechar=None)
    writer.writerow(['track name="refflat"'])

    for row in reader:
        writer.writerow(
            (
                row[2],
                row[4],
                row[5],
                row[1] + "|" + row[0],
                0,
                row[3],
                row[6],
                row[7],
                0,
                row[8],
                ",".join([str(int(x) - int(y)) for x, y in zip(row[10].split(","), row[9].split(","), strict=False) if x and y]),
                ",".join([str(int(x) - int(row[4])) for x in row[9].split(",") if x]),
            )
        )

    shutil.move(bed_output_file, bed_output_file.as_posix().lower())


def _parse_args():
    parser = argparse.ArgumentParser(
        prog="genome_generation.py",
        description="A script to generate species-specific genome-related files",
        epilog="This script uses Ensembl and UCSC provided data; this would not be possible without either of them.",
    )

    parser.add_argument("--taxon-id", type=int, dest="taxon_id")
    parser.add_argument("--release-number", type=str, dest="release_number")
    parser.add_argument("--root-save-directory", type=str, dest="root_save_directory")
    parser.add_argument("--type", type=str)
    parser.add_argument("--show-progress", action="store_true", dest="show_progress", default=False)
    return parser.parse_args()


def _main():
    args = _parse_args()

    with NCBI(taxon_id=args.taxon_id, release_number=args.release_number, fasta_type=args.type, show_progress=args.show_progress) as ncbi:
        ncbi.download_primary_assembly(save_directory=args.root_save_directory)
        ncbi.download_gtf_file(save_directory=args.root_save_directory)
        ncbi.download_cdna_fasta_file(save_directory=args.root_save_directory)
        ref_flat_file_creation(taxon_id=args.taxon_id, save_directory=args.root_save_directory)
        r_rna_interval_list_creation(taxon_id=args.taxon_id, fasta_type=args.type, save_directory=args.root_save_directory)
        bed_file_creation(taxon_id=args.taxon_id, save_directory=args.root_save_directory)


if __name__ == "__main__":
    _main()
