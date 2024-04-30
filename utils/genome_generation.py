import argparse
import csv
import ftplib
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import urllib.request
from functools import cache
from http.client import HTTPResponse

ncbi_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
ensembl_url = "ftp.ensembl.org"

ucsc_url = "https://api.genome.ucsc.edu"
ref_flat_url = "https://hgdownload.soe.ucsc.edu/goldenPath/{final_genome}/database/refFlat.txt.gz"  # fmt: skip


class Utilities:
    @cache
    @staticmethod
    def get_latest_release() -> int:
        """
        This function will identify the latest `release-###` from the ensembl_base_url

        Returns:
            int: An integer matching the release number. For example, if the latest release is `release-112`, this function returns `112`
        """
        pub_root: ftplib.FTP = ftplib.FTP(ensembl_url)
        pub_root.login()

        latest_release: int = -1
        release: str
        for release in pub_root.nlst("/pub"):
            if release.startswith("/pub/release-") and not release.endswith(".txt"):
                release_number: int = int(release.split("-")[-1])
                if release_number > latest_release:
                    latest_release = release_number

        return latest_release

    @cache
    @staticmethod
    def is_validate_release_number(release_number: str | int) -> bool:
        """
        This function will make sure `release_number` exists in the ensembl_ftp_url

        Args:
            release_number (str): The release number to check. If a string is provided, it should be in the format of (for example) `release-112`. If an integer is provided, it should only be the number (for example) `112`

        Returns:
            bool: True if the release number is valid, False otherwise
        """
        if isinstance(release_number, str):
            if release_number.startswith("release-"):
                release_number = int(release_number.split("-")[-1])
            elif release_number.isdigit():
                release_number = int(release_number)
            else:
                return False

        latest_release = Utilities.get_latest_release()
        return 19 <= release_number <= latest_release  # 19 is the earliest release number

    @cache
    @staticmethod
    def is_valid_taxon(taxon_id: int) -> bool:
        """
        This function will make sure the species provided is correct by validating it against NCBI's Taxonomy database

        Args:
            species (str): The species name

        Returns:
            bool: True if the species is valid, False otherwise
        """
        invalid_taxon_error = f"The taxonomy name you specified ({taxon_id}) is not a recognized NCBI Taxonomy name."
        response: HTTPResponse = urllib.request.urlopen(f"{ncbi_url}/taxonomy/taxon/{taxon_id}/name_report")
        as_json: dict = json.loads(response.read())

        is_valid = True
        for response_list in as_json["reports"]:
            if "errors" in response_list.keys():
                error_keys = response_list["errors"][0].keys()
                error_values = response_list["errors"][0].values()

                if ("reason" in error_keys) and (invalid_taxon_error in error_values):
                    is_valid = False

        return is_valid

    @cache
    @staticmethod
    def get_species_from_taxon(
        taxon_id: int,
        lowercase: bool = True,
        replace_spaces: bool = True,
        replace_with: str = "_",
    ) -> str:
        """
        This function will return the species name from the provided taxon ID

        It assumes the taxon ID is valid, as defined by `_is_valid_taxon`

        Args:
            taxon_id (TaxonID): The taxon ID to convert to species name
            lowercase (bool): If True, the species name will be returned in lowercase
            use_underscores (bool): If True, spaces will be replaced with underscores
            replace_with (str): The character to replace spaces with

        Returns:
            str: The species name that corresponds to the ensembl species. For example: 9606 -> homo_sapiens
        """
        response: HTTPResponse = urllib.request.urlopen(f"{ncbi_url}/taxonomy/taxon/{taxon_id}/name_report")
        as_json = json.loads(response.read())
        species = as_json["reports"][0]["taxonomy"]["current_scientific_name"]["name"]

        if lowercase:
            species = species.lower()
        if replace_spaces:
            species = species.replace(" ", replace_with)

        return species


class NCBI:
    def __init__(
        self,
        taxon_id: int,
        release_number: str,
        show_progress: bool,
    ) -> None:
        self._taxon_id: int = taxon_id
        self._release_number: str = release_number
        self._species_name: str = Utilities.get_species_from_taxon(self._taxon_id)
        self._progress: bool = show_progress

        # Validate input arguments
        if not self._release_number.startswith(("latest", "release-")):
            raise ValueError(f"""release_number must start with 'latest' or 'release-'.\nExamples: "latest", "release-112""")  # fmt: skip
        if not Utilities.is_valid_taxon(self._taxon_id):
            raise ValueError(f"The provided taxon_id ({self._taxon_id}) is not valid.\nPlease validate your taxon ID at https://www.ncbi.nlm.nih.gov/Taxonomy")  # fmt: skip
        if self._release_number != "latest" and not Utilities.is_validate_release_number(self._release_number):
            raise ValueError(f"The provided release number ({self._release_number}) is not valid. It should be in the format of `release-112` or `112`.\nPlease validate your release number at https://ftp.ensembl.org/pub")  # fmt: skip

        # Inputs are valid, we can now download the fasta primary assembly
        if self._release_number == "latest":
            self._release_number = f"release-{Utilities.get_latest_release()}"

        self._ftp = ftplib.FTP(ensembl_url)
        self._ftp.login()

    def download_fasta_file(self, save_directory: str) -> None:
        """
        This function will download the reference genome from the UCSC Genome Browser and save it as a FASTA file.

        # The full input path of the genome fasta file
        GENOME_FASTA_FILE: "genome/Mus_musculus.GRCm39.dna.primary_assembly.fa"
        GTF_FILE: "genome/Mus_musculus.GRCm39.111.gtf"
        """
        save_directory = os.path.join(save_directory, self._species_name)
        if not os.path.exists(save_directory):
            os.makedirs(save_directory, exist_ok=True)

        print(f"fasta file save directory: {save_directory}")

        total_download_size = 0
        fasta_root = f"/pub/{self._release_number}/fasta/{self._species_name}/dna"
        primary_assembly_suffix = ".dna.primary_assembly.fa.gz"
        primary_assembly_path = ""
        total_download_size: int = 0

        for filename in self._ftp.nlst(fasta_root):
            if filename.endswith(primary_assembly_suffix):
                primary_assembly_path = filename

                size = self._ftp.size(filename)
                if size is not None:
                    total_download_size += size

        # If the primary assembly path is not empty, we can download it
        if primary_assembly_path != "":
            primary_assembly_filename = primary_assembly_path.split("/")[-1]
            with open(f"{save_directory}/{primary_assembly_filename}", "wb") as fasta_file:  # fmt: skip
                self._ftp.retrbinary(f"RETR {primary_assembly_path}", fasta_file.write)
        # Otherwise we must download all *.dna.chromosome.*.fa.gz files and concatenate them
        else:
            primary_assembly_filename = ""
            chromosome_files = []
            for filename in self._ftp.nlst(fasta_root):
                if ".dna.primary_assembly." in filename:
                    chromosome_files.append(filename)
                if filename.endswith(".dna.primary_assembly.1.fa.gz"):
                    primary_assembly_filename = filename.split("/")[-1].replace(".1.fa.gz", ".fa.gz")  # fmt: skip

            # TODO: add a progress bar to this download
            with open(f"{save_directory}/{primary_assembly_filename}", "wb") as fasta_file:
                for remote_chromosome in chromosome_files:
                    chromosome_filename = remote_chromosome.split("/")[-1]
                    local_chromosome = f"{save_directory}/{chromosome_filename}"
                    with open(local_chromosome, "wb") as chr_out:
                        self._ftp.retrbinary(f"RETR {remote_chromosome}", chr_out.write)

                    with open(local_chromosome, "rb") as chr_in:
                        fasta_file.write(chr_in.read())

            # Remove primary assembly files that are not the "primary_assembly"
            for filename in os.listdir(save_directory):
                # Skip files not related to primary_assembly
                if "primary_assembly" not in filename:
                    continue
                if not filename.endswith(".dna.primary_assembly.fa.gz"):
                    os.remove(os.path.join(save_directory, filename))

        # Once the primary assembly is downloaded/created, un-gzip it
        txt_filename = primary_assembly_filename.replace(".gz", "")
        with gzip.open(f"{save_directory}/{primary_assembly_filename}", "rb") as f_in:
            with open(f"{save_directory}/{txt_filename}", "wb",) as f_out:  # fmt: skip
                shutil.copyfileobj(f_in, f_out)

        # Remove the gzipped file then convert file to lowercase name
        os.remove(f"{save_directory}/{primary_assembly_filename}")
        shutil.move(
            os.path.join(save_directory, txt_filename),
            os.path.join(save_directory, f"{self._species_name}_{self._release_number}_primary_assembly.fa"),
        )

    def download_gtf_file(self, save_directory: str) -> None:
        """
        # The full input path of the GTF file
        GTF_FILE: "genome/Mus_musculus.GRCm39.111.gtf"
        """
        save_directory = os.path.join(save_directory, self._species_name)
        if not os.path.exists(save_directory):
            os.makedirs(save_directory, exist_ok=True)

        print(f"gtf file save directory: {save_directory}")

        gtf_root = f"/pub/{self._release_number}/gtf/{self._species_name}"
        gtf_suffix = f"{self._release_number.split('-')[-1]}.gtf.gz"

        # Download file from ensembl FTP server
        # TODO: add a progress bar to this download
        gtf_filename = ""
        for filename in self._ftp.nlst(gtf_root):
            if filename.endswith(gtf_suffix):
                gtf_filename = filename.split("/")[-1]
                with open(f"{save_directory}/{gtf_filename}", "wb") as gtf_file:
                    self._ftp.retrbinary(f"RETR {filename}", gtf_file.write)

        # ungzip the gtf file
        txt_filename = gtf_filename.replace(".gz", "")
        with gzip.open(f"{save_directory}/{gtf_filename}", "rb") as f_in:
            with open(f"{save_directory}/{txt_filename}", "wb") as f_out:  # fmt: skip
                shutil.copyfileobj(f_in, f_out)

        # Remove the gzipped file then convert file to lowercase name
        os.remove(f"{save_directory}/{gtf_filename}")
        shutil.move(
            os.path.join(save_directory, txt_filename),
            os.path.join(save_directory, f"{self._species_name}_{self._release_number}.gtf"),
        )


def ref_flat_file_creation(taxon_id: int, save_directory: str, show_progress: bool) -> None:
    """
    This function will download the `refFlat.txt` file from the UCSC Genome Browser and save it as a file

    Args:
        taxon_id (TaxonID): The taxon ID to download the `refFlat.txt` file for
        save_directory (str): The directory to save the `refFlat.txt` file
    """
    species_name = Utilities.get_species_from_taxon(taxon_id)
    species_name_with_space = Utilities.get_species_from_taxon(taxon_id, replace_spaces=False)
    save_directory = os.path.join(save_directory, species_name)
    if not os.path.exists(save_directory):
        os.makedirs(save_directory, exist_ok=True)

    print(f"ref flat file save directory: {save_directory}")

    # Get a list of genomes
    response: HTTPResponse = urllib.request.urlopen(f"{ucsc_url}/list/ucscGenomes")
    as_json = json.loads(response.read())

    genome_prefix: str = ""
    genome_version: int = -1

    for genome_id in as_json["ucscGenomes"].keys():
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

    # Download `ref_flat_url` and save it as a file without creating an intermeidate gzip file
    final_genome: str = f"{genome_prefix}{genome_version}"
    download_url = ref_flat_url.format(final_genome=final_genome)
    txt_file = f"{save_directory}/ref_flat_{final_genome}.txt"

    # TODO: add a progress bar for this download
    response = urllib.request.urlopen(download_url)

    with gzip.GzipFile(fileobj=io.BytesIO(response.read()), mode="rb") as decompressed:
        with open(txt_file, "wb") as o_stream:
            o_stream.write(decompressed.read())

    shutil.move(txt_file, os.path.join(save_directory, f"{species_name}_ref_flat.txt"))


def rRNA_interval_list_creation(taxon_id: int, save_directory: str) -> None:
    """
    This function will create an interval list file suitable for CollectRnaSeqMetrics
    It should be executed at the end of the genome generation process so the required `.fa` file is availale

    Initially created by Kamil Slowikowski (December 12, 2014)
    Modified by: Arindam Ghosh (July 24, 2019)
    Modified by: Josh Loecker (April 13, 2023)

    Args:
        taxon_id (TaxonID): The taxon ID to create the rRNA interval list for
        save_directory (str): The directory to save the rRNA interval list

    Raises:
        FileNotFoundError: «The primary assembly file could not be found» if the primary assembly file is not found
    """

    species_name = Utilities.get_species_from_taxon(taxon_id)
    save_directory = os.path.join(save_directory, species_name)
    if not os.path.exists(save_directory):
        os.makedirs(save_directory, exist_ok=True)

    print(f"rRNA interval list save directory: {save_directory}")

    genome_sizes = f"{save_directory}/sizes.genome"
    genes = [file for file in os.listdir(save_directory) if file.endswith(".gtf")][0]  # fmt: skip
    primary_assembly_fa = [file for file in os.listdir(save_directory) if file.endswith(".fa")][0]  # fmt: skip
    primary_assembly_fai = f"{primary_assembly_fa}.fai"

    genes = os.path.join(save_directory, genes)
    primary_assembly_fa = os.path.join(save_directory, primary_assembly_fa)
    primary_assembly_fai = os.path.join(save_directory, primary_assembly_fai)

    rRNA_interval_list = f"{save_directory}/{species_name}_rrna.interval_list"

    # If primary_assembly_fa doesn't exist, show error and quit
    if not os.path.exists(primary_assembly_fa):
        raise FileNotFoundError(f"The primary assembly file could not be found\nSearching for: '{primary_assembly_fa}'")  # fmt: skip

    # 1. Prepare chromosome sizes file  from fasta sequence if needed.
    # Create the .fai file
    # Use a PIPE to prevent warnings from being shown
    subprocess.run(
        ["samtools", "faidx", primary_assembly_fa, primary_assembly_fai],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    primary_assembly_fai_in = open(primary_assembly_fai, "r")
    genome_sizes_out = open(genome_sizes, "w")
    for line in primary_assembly_fai_in:
        fields = line.split()
        genome_sizes_out.write(f"{fields[0]}\t{fields[1]}\n")
    primary_assembly_fai_in.close()
    genome_sizes_out.close()

    genes_in = open(genes, "r").readlines()
    genome_sizes_in = open(genome_sizes, "r").readlines()
    genome_sizes_in = [line.strip() for line in genome_sizes_in]
    rRNA_interval_list_out = open(rRNA_interval_list, "w")
    for line in genome_sizes_in:
        fields = line.split()
        rRNA_interval_list_out.write(f"@SQ\tSN:{fields[0]}\tLN:{fields[1]}\tUR:file:{primary_assembly_fa}\n")  # fmt: skip
    for line in genes_in:
        if 'gene_biotype "rRNA"' in line:
            fields = line.split()
            if fields[2] == "gene":
                gene_id = fields[-1].split('"')[1]
                rRNA_interval_list_out.write(f"{fields[0]}\t{fields[3]}\t{fields[4]}\t{fields[6]}\t{gene_id}\n")

    shutil.move(rRNA_interval_list, os.path.join(save_directory, f"{species_name}_rrna.interval_list"))

    # Sort the interval list
    # TODO: Is this needed? As of now, I don't think so
    # subprocess.run(["sort", "-k1V", "-k2n", "-k3n", "-o", rRNA_interval_list, rRNA_interval_list])  # fmt: skip


def bed_file_creation(taxon_id: int, save_directory: str) -> None:
    """
    This file will convert USCS refFlat.txt files into BED format

    # The reference BED file built from the GTF for RSEQC option
    BED_FILE: "genome/mm10_RefSeq.bed"

    modified from: informationsea at https://gist.github.com/informationsea/439d4fc53ea2b17cfb05
    """
    species_name = Utilities.get_species_from_taxon(taxon_id)
    save_directory = os.path.join(save_directory, species_name)
    bed_output_file = os.path.join(save_directory, f"{species_name}.bed")
    if not os.path.exists(save_directory):
        os.makedirs(save_directory, exist_ok=True)

    print(f"bed file save directory: {save_directory}")

    ref_flat_file = os.path.join(save_directory, f"{species_name}_ref_flat.txt")
    assert os.path.exists(ref_flat_file), f"The refFlat file could not be found\nSearching for: '{ref_flat_file}'"

    reader = csv.reader(open(ref_flat_file, "r"), delimiter="\t", quotechar=None)
    writer = csv.writer(open(bed_output_file, "w"), delimiter="\t", quotechar=None)
    writer.writerow([f'track name="refflat"'])

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
                ",".join([str(int(x) - int(y)) for x, y in zip(row[10].split(","), row[9].split(",")) if x and y]),
                ",".join([str(int(x) - int(row[4])) for x in row[9].split(",") if x]),
            )
        )

    shutil.move(bed_output_file, bed_output_file.lower())


def parse_args():
    parser = argparse.ArgumentParser(
        prog="genome_generation.py",
        description="A script to generate species-specific genome-related files",
        epilog="This script uses Ensembl and UCSC provided data; this would not be possible without either of them.",
    )

    parser.add_argument("--taxon-id", type=int, dest="taxon_id")
    parser.add_argument("--release-number", type=str, dest="release_number")
    parser.add_argument("--root-save-directory", type=str, dest="root_save_directory")
    parser.add_argument("--show-progress", action="store_true", dest="show_progress", default=False)

    args = parser.parse_args()
    return args


def main():
    # fmt: off
    args = parse_args()

    ncbi = NCBI(taxon_id=args.taxon_id, release_number=args.release_number, show_progress=args.show_progress)
    ncbi.download_fasta_file(save_directory=args.root_save_directory)
    ncbi.download_gtf_file(save_directory=args.root_save_directory)

    ref_flat_file_creation(taxon_id=args.taxon_id, save_directory=args.root_save_directory, show_progress=args.show_progress)
    rRNA_interval_list_creation(taxon_id=args.taxon_id, save_directory=args.root_save_directory)
    bed_file_creation(taxon_id=args.taxon_id, save_directory=args.root_save_directory)


if __name__ == "__main__":
    main()
