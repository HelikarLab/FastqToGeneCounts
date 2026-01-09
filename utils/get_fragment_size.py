"""Calculate fragment size for each gene/transcript.

For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

import re
import sys
from collections.abc import Sequence
from optparse import OptionParser
from pathlib import Path
from statistics import mean, median, pstdev

import pysam
from tqdm import tqdm

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "5.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def overlap_length2(exons: Sequence[Sequence[int]], read_start: int, next_ref_start: int) -> int:
    """Sum the overlap between sorted, non-overlapping exon intervals with inclusive coordinates.

    :param exons: Sorted by start, non-overlapping: [(start_inclusive, end_inclusive), ...]
    :param read_start: The start of the potential overlapping read
    :param next_ref_start: The start if the read's paired mate
    :return: The total overlap length
    """
    total = 0
    for x1, x2 in exons:
        if x2 <= read_start:  # exon is entirely before [start, end); skip it
            continue

        # exon starts at/after the next reference start
        # because the exons are sorted, no exons after this can overlap
        if x1 >= next_ref_start:
            break

        # compute the overlap of [x1, x2) with [a, b)
        lo = max(read_start, x1)
        hi = min(next_ref_start, x2)
        if lo < hi:
            total += hi - lo

    return total


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


def _fragment_size(reference_bed_filepath: Path, bam_filepath: Path, qcut: int, ncut: int, threads: int):
    """Calculate the fragment size for each gene."""
    ref_size_bytes = reference_bed_filepath.stat().st_size
    with (
        reference_bed_filepath.open("r") as bed_stream,
        pysam.AlignmentFile(bam_filepath.as_posix(), threads=threads) as alignment_file,
        tqdm(total=ref_size_bytes, desc=f"Calculating fragment sizes: '{bam_filepath.name}'", unit="B", unit_scale=True, unit_divisor=1024) as pbar,
    ):
        all_sizes = []
        references = set(alignment_file.references)
        for i, line in enumerate(bed_stream):
            if i == 100:
                print(all_sizes)
                break
            line_bytes = len(line)
            if line.startswith(("#", "track", "browser")):
                pbar.update(line_bytes)
                continue
            fields = line.split()
            chrom = fields[0]
            if chrom not in references:
                alt_chrom = chrom.removeprefix("chr") if chrom.startswith("chr") else f"chr{chrom}"
                if alt_chrom in references:
                    chrom = alt_chrom
                else:
                    pbar.update(line_bytes)
                    continue  # skip this record; contig doesn't exist in BAM file

            tx_start = int(fields[1])
            tx_end = int(fields[2])
            gene_name = fields[3]

            exon_starts: list[int] = _calculate_exon_start(tx_start, fields)
            exon_ends: list[int] = _calculate_exon_end(exon_starts, fields)
            gene_id: str = "\t".join([str(i) for i in (chrom, tx_start, tx_end, gene_name)])
            exon_range: list[tuple[int, int]] = list(zip(exon_starts, exon_ends, strict=True))

            frag_sizes = []
            aligned_reads = alignment_file.fetch(contig=chrom, start=tx_start, stop=tx_end)
            for read in aligned_reads:
                if (
                    not read.is_paired  # single-end sequencing, we only want paired-end
                    or read.is_read2  # reverse read, we want forward read
                    or read.mate_is_unmapped  # paired read is not mapped
                    or read.is_qcfail  # low quality
                    or read.is_duplicate  # duplicate read
                    or read.is_secondary  # non primary hit
                    or read.mapping_quality < qcut  # low mapping quality
                ):
                    continue

                read_start = read.reference_start
                next_ref_start = read.next_reference_start
                if read_start > next_ref_start:
                    read_start, next_ref_start = next_ref_start, read_start
                if read_start < tx_start or next_ref_start > tx_end:
                    continue
                
                length = overlap_length2(exon_range, read_start=read_start, next_ref_start=next_ref_start) + read.query_alignment_length
                frag_sizes.append(length)
            
            all_sizes.extend(frag_sizes)

            yield (
                "\t".join([str(i) for i in (gene_id, len(frag_sizes), 0, 0, 0)])
                if len(frag_sizes) < ncut
                else "\t".join([str(i) for i in (gene_id, len(frag_sizes), mean(frag_sizes), median(frag_sizes), pstdev(frag_sizes))])
            )
            pbar.update(line_bytes)


def main():
    usage = "%prog [options]" + "\n" + (__doc__ or "") + "\n"
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
        for tmp in _fragment_size(
            reference_bed_filepath=options.refgene_bed,
            bam_filepath=options.input_file,
            qcut=options.map_qual,
            ncut=options.fragment_num,
            threads=4,
        ):
            o_stream.write(tmp + "\n")


if __name__ == "__main__":
    main()
