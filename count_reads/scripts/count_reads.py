#!/usr/bin/env python3
"""
Generate count tables from sequence data.

Read (paired) fasta/fastq files and count the occurance of each read or sub-region thereof.
Gzip compression is supported and filetype is detected via the file extension.
Observed reads can be compared to an expected library, including fuzzy matching.
Outputs a TSV count table.

Regions are defined in a TSV with four columns. Lines starting with # are skipped:
Name - The region name
Read - Match the reverse (entry is any leading substring of reverse) or forward (otherwise) read
Start - 0-based inclusive start of the region
End - 0-based exclusive end of the read

The library to compare to is also defined in a TSV file, with one column per region to be compared to the library. Each row contains the sequences of these regions from one oligo (the combination observed in reads will also be checked). Optional per region max accepted distances can be given at the end of the headers after a bar (|).

Output can be grouped by features of the sequence name by supplying a RegEx with a capture group
to the --group argument. The result of the capture group is used to group output counts as if it
were an additional region.
"""
import argparse
import warnings
import sys
import multiprocessing
import re
from itertools import islice

from tqdm import tqdm
from count_reads import *

def main():
    """
    Main
    """
    # from bin.count_reads import *
    # args = parse_args(["--library", "config/swallow.library", "--regions", "config/single_cell_swallow.regions", "--min_phred", "10", "--sort", "--threads", "1", "--verbose", "data/scdna_seq1/hap1_54-2-1/chrl.trimmed.fastq"])
    args = parse_args()

    # Setup progress bar style (slower \n on redirected stderr, doesn't detect ipython)
    if not args.verbose:
        bar = create_default_tqdm(disable=True)
    elif not sys.stderr.isatty():
        bar = create_default_tqdm(newline=True, mininterval=5, maxinterval=60)
    else:
        bar = tqdm

    regions = Region.parse_regions(args.regions, progress=bar) if args.regions else None
    library = SequenceLibrary.from_file(args.library, progress=bar) if args.library else None

    group = re.compile(args.group) if args.group else None

    seqs = SequenceParser(args.file1, args.file2, filetype=args.format,
                          gzipped=args.gzip, reverse_complement=True)

    # Only run up to 5k sequences in debug mode
    if args.debug:
        seqs = islice(seqs, 5000)

    # reset progress bar for writing to log to avoid many similar lines
    if args.verbose and not sys.stderr.isatty():
        bar = create_default_tqdm(newline=True, mininterval=60, maxinterval=300)

    counts, filtered_seqs = count_seqs(seqs, regions, group=group,
                                       phred_threshold=args.min_phred,
                                       progress=bar)
    if filtered_seqs.total > 0:
        print(str(filtered_seqs), file=sys.stderr)

    if library is not None:
        # adjust progress bar for writing to log
        if args.verbose and not sys.stderr.isatty():
            bar = create_default_tqdm(newline=True, mininterval=10, miniters=int(0.01*len(counts)))
        counts = compare_to_library(counts, library, threshold=args.max_distance,
                                    threads=args.threads, progress=bar, chunksize=args.chunksize)

    write_counts(counts, sort=args.sort, progress=bar if not sys.stdout.isatty() else False)

def arg_parser():
    """
    Construct argument parser
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("file1", metavar="F1", help="File to count sequences from")

    parser.add_argument("file2", metavar="F2", nargs="?", help="File to count sequences from")

    parser.add_argument("--regions", "-r", help="Path to TSV file defining regions of interest")

    parser.add_argument("--library", "-l", help="Path to TSV file defining expected library")

    parser.add_argument("--group", "-g", help="Group output by this feature")

    parser.add_argument("--sort", "-s", action="store_true", help="Sort output by descending count")

    parser.add_argument("--min_phred", "-p", type=int, default=0,
                        help="Filter reads with mean Phred score below this threshold")

    parser.add_argument("--max_distance", "-d", type=int, default=3,
                        help="Only assign nearest sequences with at most this many substitutions")

    parser.add_argument("--format", "-f", choices=["auto", "fasta", "fastq"],
                        default="auto", help="Input file(s) format")

    parser.add_argument("--gzip", "-z", choices=["auto", "gzip", "none"],
                        default="auto", help="Input file(s) is gzipped")

    parser.add_argument("--threads", "-t", default=1, type=int, help="Number of threads to use")

    parser.add_argument("--chunksize", "-c", type=int, help="Chunksize for parallel processing")

    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

    parser.add_argument("--debug", action="store_true", help="Debug mode")

    return parser

def parse_args(arg_list=None):
    """
    Parse and validate script arguments
    """
    args = arg_parser().parse_args(arg_list)

    if args.min_phred < 0 or args.min_phred > 42:
        warnings.warn(("Standard Phred scores range from 0 (low qaulity) to 42 (high quality). "
                       f"{args.min_phred} is outside this range."))

    if args.threads > multiprocessing.cpu_count():
        warnings.warn((f"More threads allocated ({args.threads}) than visible CPUs "
                       f"({multiprocessing.cpu_count()})"))

    if args.chunksize is not None and args.threads == 1:
        warnings.warn("--chunksize set but --threads == 1, will have no effect")

    return args

if __name__ == "__main__":
    main()
