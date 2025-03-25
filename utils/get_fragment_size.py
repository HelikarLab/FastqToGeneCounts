"""This script is taken from: https://rseqc.sourceforge.net/.

Specifically, the `scripts/RNA_fragment_size.py` file

calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

import os
import sys
from optparse import OptionParser

import pysam
from numpy import mean, median, std

if sys.version_info[0] != 3:
    print(
        "\nYou are using python" + str(sys.version_info[0]) + "." + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n",
        file=sys.stderr,
    )
    sys.exit()


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "5.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def overlap_length2(lst1, lst2):
    l = 0
    for x in lst1:
        for y in lst2:
            l += len(list(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)))
    return l


def fragment_size(bedfile, samfile, qcut=30, ncut=5):
    """Calculate the fragment size for each gene"""
    for line in open(bedfile):
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

        for st, end in zip(exon_starts, exon_ends, strict=False):
            exon_range.append([st + 1, end + 1])
        # exon_range.append([chrom, st,end])

        try:
            alignedReads = samfile.fetch(chrom, tx_start, tx_end)
        except:
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
            yield "\t".join([str(i) for i in (geneID, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))])


def main():
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

    if not os.path.exists(options.input_file):
        raise FileNotFoundError(f"BAM file does not exist: '{options.input_file}'")
    if not os.path.exists(options.bai_file):
        raise FileNotFoundError(f"BAI file does not exist: '{options.bai_file}'")
    if not os.path.exists(options.refgene_bed):
        raise FileNotFoundError(f"BED file does not exist: '{options.refgene_bed}'")

    os.makedirs(os.path.basename(options.output_filepath), exist_ok=True)
    with open(options.output_filepath, "w") as o_stream:
        o_stream.write("\t".join([str(i) for i in ("chrom", "tx_start", "tx_end", "symbol", "frag_count", "frag_mean", "frag_median", "frag_std")]))
        o_stream.write("\n")
        for tmp in fragment_size(options.refgene_bed, pysam.Samfile(options.input_file), options.map_qual, options.fragment_num):
            o_stream.write(tmp + "\n")


if __name__ == "__main__":
    main()
