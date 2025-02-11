"""
This script is taken from: https://rseqc.sourceforge.net/

calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

import sys
from typing import TextIO
from pathlib import Path
import statistics

from snakemake.script import snakemake
import pysam


if sys.version_info[0] != 3:
    print(
        "\nYou are using python"
        + str(sys.version_info[0])
        + "."
        + str(sys.version_info[1])
        + " This verion of RSeQC needs python3!\n",
        file=sys.stderr,
    )
    sys.exit()

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def overlap_length2(lst1, lst2):
    length = 0
    for x in lst1:
        for y in lst2:
            length += len(list(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)))
    return length


def fragment_size(bed_filepath: Path, sam_filepath: Path, qcut=30, ncut=5):
    """calculate the fragment size for each gene"""

    sam_obj = pysam.AlignmentFile(sam_filepath)
    for line in bed_filepath.open(mode="r"):
        if line.startswith(("#", "track", "browser")):
            continue

        fields = line.split()

        chrom = fields[0]
        tx_start = int(fields[1])
        tx_end = int(fields[2])
        gene_name = fields[3]
        strand = fields[5].replace(" ", "_")

        exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
        exon_starts = list(map((lambda x: x + tx_start), exon_starts))
        exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
        exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
        gene_id = "\t".join([str(i) for i in (chrom, tx_start, tx_end, gene_name)])
        exon_range = [[start + 1, end + 1] for start, end in zip(exon_starts, exon_ends)]

        try:
            aligned_reads = sam_obj.fetch(chrom, tx_start, tx_end)
        except ValueError:
            yield "\t".join([str(i) for i in (gene_id, 0, 0, 0)])
            continue

        frag_sizes = []
        for aligned_read in aligned_reads:
            if (
                not aligned_read.is_paired  # skip single sequencing
                or aligned_read.is_read2
                or aligned_read.mate_is_unmapped
                or aligned_read.is_qcfail  # skip low quality
                or aligned_read.is_duplicate  # skip duplicate read
                or aligned_read.is_secondary  # skip non-primary hit
                or aligned_read.mapping_quality < qcut
            ):
                continue

            read_start = aligned_read.reference_start
            mate_start = aligned_read.next_reference_start
            if read_start > mate_start:
                (read_start, mate_start) = (mate_start, read_start)
            if read_start < tx_start or mate_start > tx_end:
                continue
            read_len = aligned_read.query_length
            map_range = [[read_start + 1, mate_start]]
            frag_len = overlap_length2(exon_range, map_range) + read_len
            frag_sizes.append(frag_len)

        mean = statistics.mean(frag_sizes)
        median = statistics.median(frag_sizes)
        std = statistics.stdev(frag_sizes)
        yield (
            "\t".join([str(i) for i in (gene_id, len(frag_sizes), 0, 0, 0)])
            if len(frag_sizes) < ncut
            else "\t".join([str(i) for i in (gene_id, len(frag_sizes), mean, median, std)])
        )


def write_header(o_stream: TextIO):
    o_stream.write(
        "\t".join(
            [
                str(i)
                for i in (
                    "chrom",
                    "tx_start",
                    "tx_end",
                    "symbol",
                    "frag_count",
                    "frag_mean",
                    "frag_median",
                    "frag_std",
                )
            ]
        )
    )
    o_stream.write("\n")


def main():
    print("\n")

    bam_file_input = Path(snakemake.input[0])
    bai_file_input = Path(snakemake.input[1])
    bed_file_input = Path(snakemake.input[2])
    output_filepath = Path(snakemake.output[0])

    if not bam_file_input.exists():
        raise FileNotFoundError(f"Input file does not exist: '{bam_file_input}'")
    if not bai_file_input.exists():
        raise FileNotFoundError(f"Input file does not exist: '{bai_file_input}'")
    if not bed_file_input.exists():
        raise FileNotFoundError(f"Input file does not exist: '{bed_file_input}'")

    if not bam_file_input.suffix == ".bam":
        raise ValueError(f"Input file must be a BAM file, got file: '{bam_file_input}'")
    if not bai_file_input.suffix == ".bai":
        raise ValueError(f"Input file must be a BAI file, got file: '{bai_file_input}'")
    if not bed_file_input.suffix == ".bed":
        raise ValueError(f"Input file must be a BED file, got file: '{bed_file_input}'")

    num_lines = sum(1 for _ in bed_file_input.read_text())
    with output_filepath.open(mode="w") as o_stream:
        write_header(o_stream)
        for i, line in enumerate(
            fragment_size(
                bed_filepath=bed_file_input,
                sam_filepath=bam_file_input,
                qcut=30,
                ncut=3,
            )
        ):
            # Print progress every 10%
            if i % (num_lines // 10) == 0:
                print(f"Progress: {i / num_lines * 100:.0f}%")

            o_stream.write(f"{line}\n")

    print("\n")


if __name__ == "__main__":
    main()
