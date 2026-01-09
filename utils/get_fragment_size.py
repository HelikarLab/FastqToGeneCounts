"""Calculate fragment size for each gene/transcript.

For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

import re
import sys
from optparse import OptionParser
from pathlib import Path

import pysam
from numpy import mean, median, std

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "5.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def overlap_length2(lst1, lst2):
    """Calculate the overlap length between two lists of intervals."""
    length = 0
    for x in lst1:
        for y in lst2:
            length += len(list(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)))
    return length


def _calculate_exon_start(tx_start: int, fields: list[str]) -> list[int]:
    exon_starts: list[int] = list(map(int, fields[11].rstrip(",\n").split(",")))
    return [i + tx_start for i in exon_starts]


def _calculate_exon_end(exon_starts: list[int], fields: list[str]) -> list[int]:
    exon_ends: list[int] = list(map(int, fields[10].rstrip(",\n").split(",")))
    return [x + y for x, y in zip(exon_starts, exon_ends, strict=True)]


def get_contig(chrom_value: str) -> str:
    chromosome_match = re.match(r"^chr(\d+)$", chrom_value)
    if chromosome_match:
        return chromosome_match.group(1)
    if chrom_value.endswith(("X", "Y")):
        return chrom_value[-1]

    if "_" in chrom_value:
        split_chrom = chrom_value.split("_")
        contig = split_chrom[-2] if split_chrom[-1] == "random" else split_chrom[-1]
        contig = contig.replace("v", ".")
        return contig

    raise ValueError(f"Cannot parse chromosome name: {chrom_value}")


def fragment_size(reference_bed_filepath: Path, input_bam_filepath: Path, qcut: int = 30, ncut: int = 5):
    """Calculate the fragment size for each gene."""
    with reference_bed_filepath.open("r") as i_stream, pysam.AlignmentFile(input_bam_filepath.as_posix()) as sam_file:
        for line in i_stream:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            chrom = fields[0]
            contig = get_contig(chrom_value=chrom)

            tx_start = int(fields[1])
            tx_end = int(fields[2])
            gene_name = fields[3]

            exon_starts: list[int] = _calculate_exon_start(tx_start, fields)
            exon_ends = _calculate_exon_end(exon_starts, fields)
            gene_id = "\t".join([str(i) for i in (chrom, tx_start, tx_end, gene_name)])
            exon_range = [[st + 1, end + 1] for st, end in zip(exon_starts, exon_ends, strict=True)]

            frag_sizes = []
            aligned_reads = sam_file.fetch(contig=get_contig(chrom_value=chrom), start=tx_start, stop=tx_end)
            for read in aligned_reads:
                if (
                    not read.is_paired  # single-end sequencing
                    or read.is_read2  # reverse read
                    or read.mate_is_unmapped  # paired read is not mapped
                    or read.is_qcfail  # low quality
                    or read.is_duplicate  # duplicate read
                    or read.is_secondary  # non primary hit
                    or read.mapping_quality < qcut  # low mapping quality
                ):
                    continue

                read_start = read.reference_start
                mate_start = read.next_reference_start
                if read_start > mate_start:
                    read_start, mate_start = mate_start, read_start
                if read_start < tx_start or mate_start > tx_end:
                    continue
                map_range = [[read_start + 1, mate_start]]
                frag_sizes.append(overlap_length2(exon_range, map_range) + read.query_alignment_length)
            yield (
                "\t".join([str(i) for i in (gene_id, len(frag_sizes), 0, 0, 0)])
                if len(frag_sizes) < ncut
                else "\t".join([str(i) for i in (gene_id, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))])
            )


def _main():
    usage = "%prog [options]" + "\n" + __doc__ + "\n"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-i", "--input", action="store", type="string", dest="input_file", help="Input BAM file")
    parser.add_option(
        "-r",
        "--refgene",
        action="store",
        type="string",
        dest="refgene_bed",
        help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]",
    )
    parser.add_option("-a", "--bai", action="store", type="string", dest="bai_file", help="Input BAI file that matches input BAM file")
    parser.add_option(
        "-q",
        "--mapq",
        action="store",
        type="int",
        dest="map_qual",
        default=30,
        help='Minimum mapping quality (phred scaled) for an alignment to be called "uniquely mapped". default=%default',
    )
    parser.add_option(
        "-n", "--frag-num", action="store", type="int", dest="fragment_num", default=3, help="Minimum number of fragment. default=%default"
    )
    parser.add_option("-o", "--output", action="store", type="string", dest="output_filepath", help="The output filepath.")

    (options, args) = parser.parse_args()

    if not (options.input_file and options.refgene_bed):
        parser.print_help()
        sys.exit(0)

    options.input_file = Path(options.input_file)
    options.refgene_bed = Path(options.refgene_bed)
    options.bai_file = Path(options.bai_file)
    options.output_filepath = Path(options.output_filepath)

    if not options.input_file.exists():
        raise FileNotFoundError(f"BAM file does not exist: '{options.input_file}'")
    if not options.bai_file.exists():
        raise FileNotFoundError(f"BAI file does not exist: '{options.bai_file}'")
    if not options.refgene_bed.exists():
        raise FileNotFoundError(f"BED file does not exist: '{options.refgene_bed}'")

    with options.output_filepath.open("w") as o_stream:
        o_stream.write("\t".join([str(i) for i in ("chrom", "tx_start", "tx_end", "symbol", "frag_count", "frag_mean", "frag_median", "frag_std")]))
        o_stream.write("\n")
        for tmp in fragment_size(
            reference_bed_filepath=options.refgene_bed,
            input_bam_filepath=options.input_file,
            qcut=options.map_qual,
            ncut=options.fragment_num,
        ):
            o_stream.write(tmp + "\n")


if __name__ == "__main__":
    _main()
