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

import snakemake
from numpy import mean, median, std

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


def fragment_size(bedfile, samfile, qcut=30, ncut=5):
    """calculate the fragment size for each gene"""
    for line in open(bedfile, "r"):
        exon_range = []
        if line.startswith(("#", "track", "browser")):
            continue
        fields = line.split()

        chrom = fields[0]
        tx_start = int(fields[1])
        tx_end = int(fields[2])
        geneName = fields[3]
        trand = fields[5].replace(" ", "_")

        exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
        exon_starts = list(map((lambda x: x + tx_start), exon_starts))
        exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
        exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
        geneID = "\t".join([str(i) for i in (chrom, tx_start, tx_end, geneName)])

        for st, end in zip(exon_starts, exon_ends):
            exon_range.append([st + 1, end + 1])

        try:
            alignedReads = samfile.fetch(chrom, tx_start, tx_end)
        except:  # noqa: E722
            yield "\t".join([str(i) for i in (geneID, 0, 0, 0)])
            continue

        frag_sizes = []
        for aligned_read in alignedReads:
            if not aligned_read.is_paired:  # skip single sequencing
                continue
            if aligned_read.is_read2:
                continue
            if aligned_read.mate_is_unmapped:
                continue
            if aligned_read.is_qcfail:
                continue  # skip low quanlity
            if aligned_read.is_duplicate:
                continue  # skip duplicate read
            if aligned_read.is_secondary:
                continue  # skip non primary hit
            if aligned_read.mapq < qcut:
                continue

            read_st = aligned_read.pos
            mate_st = aligned_read.pnext
            if read_st > mate_st:
                (read_st, mate_st) = (mate_st, read_st)
            if read_st < tx_start or mate_st > tx_end:
                continue
            read_len = aligned_read.qlen
            map_range = [[read_st + 1, mate_st]]
            # map_range = [[chrom, read_st, mate_st]]
            frag_len = overlap_length2(exon_range, map_range) + read_len
            frag_sizes.append(frag_len)
        if len(frag_sizes) < ncut:
            yield "\t".join([str(i) for i in (geneID, len(frag_sizes), 0, 0, 0)])
        else:
            yield "\t".join(
                [str(i) for i in (geneID, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))]
            )


def write_header(o_stream: TextIO):
    o_stream.write(
        "\t".join(
            [
                str(i)
                for i in ("chrom", "tx_start", "tx_end", "symbol", "frag_count", "frag_mean", "frag_median", "frag_std")
            ]
        )
    )
    o_stream.write("\n")


def main():
    print("\n")
    num_lines = sum(1 for i in open(snakemake.config["BED_FILE"]))  # noqa
    with open(snakemake.output[0], "w") as o_stream:  # noqa
        write_header(o_stream)
        for i, line in enumerate(
            fragment_size(bedfile=snakemake.config["BED_FILE"], samfile=snakemake.input.bam, qcut=30, ncut=3)  # noqa
        ):
            # Print progress every 10%
            if i % (num_lines // 10) == 0:
                print(f"Progress: {i / num_lines * 100:.0f}%")

            o_stream.write(f"{line}\n")

    print("\n")


if __name__ == "__main__":
    main()
