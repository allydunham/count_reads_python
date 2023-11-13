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
import pathlib
import gzip
import sys
import multiprocessing
import re
from functools import partial
from itertools import repeat, product, islice
from dataclasses import dataclass

import numpy as np
from Bio import SeqIO, SeqRecord
from tqdm import tqdm
from tqdm.utils import _unicode, disp_len

# TODO tests
# TODO optimise
# TODO collapse results to only nearest
# TODO Only library results
# TODO Only unique matches
# TODO summary statistics (recombination etc.)
# TODO Translate into stand-alone tool (probably in rust/c++, maybe with python bindings?)

FA_EXTS = {".fa": "fasta", ".fasta": "fasta", ".fastq": "fastq", ".fq": "fastq"}

@dataclass
class ReadPair:
    """Pair of sequencing reads, with a forward and reverse Bio.SeqRecord"""
    forward: SeqRecord.SeqRecord
    reverse: SeqRecord.SeqRecord = None

    @property
    def phred_scores(self):
        return (self.forward.letter_annotations['phred_quality'],
                self.reverse.letter_annotations['phred_quality'] if self.reverse else np.nan)

class SequenceParser:
    """Generator returning readpairs from a fasta or pair of fasta files"""
    def __init__(self, path1, path2=None, filetype="auto", gzipped="auto", reverse_complement=True):
        self.path1 = path1
        self.path2 = path2
        self.filetype = filetype
        self.gzipped = gzipped
        self.reverse_complement = reverse_complement
        self.paired = path2 is not None
        self._len = None

    def __repr__(self):
        return (f"SequenceParser({self.path1!r}, {self.path2!r}, {self.filetype!r}, "
                f"{self.gzipped!r}, {self.reverse_complement!r})")

    def __iter__(self):
        fasta1 = self.parse_fasta(self.path1, self.filetype, self.gzipped)
        fasta2 = self.parse_fasta(self.path2, self.filetype, self.gzipped) if self.path2 is not None else repeat(None)

        for seq1, seq2 in zip(fasta1, fasta2):
            if seq2 is not None and self.reverse_complement:
                seq2 = seq2.reverse_complement(id=True, name=True, description=True)
            yield ReadPair(seq1, seq2)

    def __len__(self):
        # Cache length neively - this can go wrong if the file changes between calls, but that isn't expected behaviour
        if self._len is None:
            self._len = sum(1 for _ in self)
        return self._len

    @staticmethod
    def parse_fasta(path, filetype="auto", gzipped="auto"):
        """
        Generator yielding sequences from a Fasta/Fastq file, handling file type and gzip compression
        """
        exts = pathlib.Path(path).suffixes
        if gzipped == "auto" and exts[-1] == ".gz" or gzipped == "gzip":
            _open = partial(gzip.open, mode="rt")
            exts = exts[:-1]
        else:
            _open = open

        mode = FA_EXTS[exts[-1]] if filetype == "auto" else filetype

        with _open(path) as file:
            for seq in SeqIO.parse(file, mode):
                yield seq

class SequenceLibrary:
    """Library of expected sequence regions"""
    def __init__(self, max_distances=None, **kwargs):
        self.regions = {k: set(v) for k, v in kwargs.items()}
        self.max_distances = (
            {k: None for k in self.regions} if max_distances is None else max_distances
        )

        if self.max_distances.keys() != self.regions.keys():
            raise ValueError("Distance keys do not match region names")

        self.matrices = {
            k: np.array([[base for base in seq] for seq in self.regions[k]]) for k in kwargs.keys()
        }
        self.combinations = set(i for i in zip(*kwargs.values()))

    def __repr__(self):
        seqs = zip(*self.combinations)
        args = ", ".join(f"{name}={seqs!r}" for name, seqs in zip(self.regions.keys(), seqs))
        return f"SequenceLibrary({args})"

    def __getitem__(self, x):
        return self.regions[x]

    @property
    def keys(self):
        return self.regions.keys()

    @classmethod
    def from_file(cls, path, progress=False):
        """Read a sequence library from a TSV file"""
        progress = _get_progress(progress)
        library = {}
        with open(path) as lib_file:
            headers = next(lib_file).strip().split()

            # Extract distances from headers
            distances = {}
            for i, h in enumerate(headers):
                if "|" in h:
                    h = h.split("|")
                    distances[h[0]] = int(h[1])
                    headers[i] = h[0]
                else:
                    distances[h] = None

            for h in headers:
                library[h] = []

            for line in progress(lib_file, "Parsing library"):
                line = line.strip().upper().split()
                for i,s in enumerate(line):
                    library[headers[i]].append(s)

        return cls(**library)

class Region:
    """
    Region of sequence read
    """
    def __init__(self, name, start, end, reverse=False):
        self.name = name
        self.start = start
        self.end = end
        self.reverse = reverse

    def __repr__(self):
        return f"Region({self.name!r}, {self.start!r}, {self.end!r}, reverse={self.reverse!r})"

    def __call__(self, readpair):
        seq = readpair.reverse if self.reverse else readpair.forward
        if seq is None:
            raise ValueError((f"{'Reverse' if self.reverse else 'Forward'} read missing "
                              "from read pair"))

        start = len(seq) + self.start if self.start < 0 else self.start
        end = len(seq) if self.end is None else self.end

        seq_slice = str(seq[start:end].seq)

        if len(seq_slice) != end - start:
            return None
        return ObservedRegion(self.name, seq_slice, self.start, self.end, self.reverse)

    @classmethod
    def permissive_region(cls, reverse=False):
        """
        Return a region permitting the entire sequence
        """
        return cls("reverse" if reverse else "forward", 0, None, reverse)

    @classmethod
    def from_line(cls, line):
        """
        Generate a Region from a TSV configuration line
        """
        line = line.split("\t")
        reverse =  "reverse".startswith(line[1].lower()) # Partial match reverse
        end = None if line[3].lower() in ("end", "none") else int(line[3])
        return cls(line[0], int(line[2]), end, reverse=reverse)

    @classmethod
    def parse_regions(cls, path, progress=False):
        """
        Parse regions from a config TSV file, returning a list of Region objects
        """
        progress = _get_progress(progress)

        regions = []
        with open(path) as region_file:
            for line in progress(region_file, desc="Reading regions"):
                if line[0] == "#":
                    continue
                regions.append(cls.from_line(line.strip()))

        names = [r.name for r in regions]
        if len(set(names)) < len(names):
            raise ValueError("Duplicate region names - regions must be uniquely named")

        return regions

class ObservedRegion:
    """Observed subsection of a read with library information"""
    def __init__(self, name, sequence, start=None, end=None, reverse=None,
                 nearest=None, distance=None, in_library=None):
        self.name = name
        self.sequence = sequence
        self.start = start
        self.end = end
        self.reverse = reverse
        self.nearest = nearest
        self.distance = distance
        self.in_library = in_library

    def __repr__(self):
        return (f"ObservedRegion({self.name!r}, {self.sequence!r}, {self.start!r}, {self.end!r}, "
                f"{self.reverse!r}, {self.nearest!r}, {self.distance!r}, {self.in_library!r})")

    def compare_library(self, library, threshold=3):
        """
        Determine the nearest sequence from a known library

        library: numpy array of library sequences with one column per bp
        threshold: maximum hamming distance to consider nearby
        """
        if not len(self.sequence) == library.shape[1]:
            raise ValueError((f"Library sequences (lenth {library.shape[1]}) don't match "
                              f"sequence (length {len(self.sequence)})"))

        # Calculate nearest sequence by hamming distance
        seq = np.array([i for i in self.sequence])

        hamming = np.count_nonzero(library != seq, axis=1) # Levenstein better?
        min_hamming = np.min(hamming)

        self.distance = int(min_hamming)
        self.in_library = self.distance == 0

        if min_hamming <= threshold:
            self.nearest = [
                "".join(i) for i in library[np.argwhere(hamming == min_hamming).flatten()]
            ]

class ObservedCombination:
    """
    Container for observed combination of region sequences with read count and library information
    """
    def __init__(self, regions, count=0, group=None, in_library=None, nearest_in_library=None):
        self.regions = regions
        self.count = count
        self.group = group
        self.in_library = in_library
        self.nearest_in_library = nearest_in_library

    def __repr__(self):
        return (f"ObservedCombination({self.regions!r}, {self.count!r}, {self.group!r},"
                f"{self.in_library!r}, {self.nearest_in_library!r})")

    @staticmethod
    def get_key(regions, group=None):
        """Generate a key for paired regions and sequences"""
        regs = "-".join(f"{r.name}:{r.sequence}" for r in regions.values())
        if group is None:
            return regs
        return f"group:{group}-{regs}"

    def compare_library(self, library, threshold=3):
        """Determine the nearest library sequence for each observed region"""
        for k, r in self.regions.items():
            if k in library.keys:
                t = library.max_distances[k] if library.max_distances[k] is not None else threshold
                r.compare_library(library.matrices[k], threshold=t)

        comb_seq = tuple(self.regions[i].sequence for i in library.keys)
        self.in_library = comb_seq in library.combinations
        nearest_seqs = [self.regions[i].nearest for i in library.keys]
        if all(i is not None for i in nearest_seqs):
            self.nearest_in_library = sum(i in library.combinations for i in product(*nearest_seqs))

class FilterCount:
    def __init__(self):
        self.total = 0
        self.forward = 0
        self.reverse = 0
        self.length = 0

    def __str__(self):
        return (f"{self.total} filtered reads\n"
                f"\t{self.forward} low quality forward reads\n"
                f"\t{self.reverse} low quality reverse reads\n"
                f"\t{self.length} incorrect length reads")

    def add_error(self, forward=False, reverse=False, length=False):
        """Add error counts"""
        self.total += 1
        self.forward += forward
        self.reverse += reverse
        self.length += length

class _MultithreadCompareWrapper:
    """Wrapper to call .compare_library on a ObservedCombination with a library/threshold"""
    def __init__(self, library, threshold):
        self.library = library
        self.threshold = threshold

    def __call__(self, x):
        x.compare_library(self.library, threshold=self.threshold)
        return x

def count_seqs(seqs, regions=None, group=None, phred_threshold=0, progress=None):
    """
    Count occurance of read forms with option to compare each to a library

    Returns a list of ObservedCombinations and a FilterCount of filtered reads
    """
    progress = _get_progress(progress)
    filtered_seqs = FilterCount()
    counts = {}

    if regions is None:
        regions = [Region.permissive_region(reverse=False)]
        if seqs.paired:
            regions.append(Region.permissive_region(reverse=True))

    for readpair in progress(seqs, total=float("inf"), desc="Counting sequence combinations"):
        # Check read quality
        phreds = [np.mean(x) for x in readpair.phred_scores]

        if phreds[0] < phred_threshold or phreds[1] < phred_threshold:
            filtered_seqs.add_error(forward = phreds[0] < phred_threshold,
                                    reverse = phreds[1] < phred_threshold)
            continue

        if group is not None:
            read_group = group.search(readpair.forward.description)
            read_group = read_group.group(1) if read_group else ""
        else:
            read_group = None

        seq_regions = {r.name: r(readpair) for r in regions}

        if any(i is None for i in seq_regions.values()):
            filtered_seqs.add_error(length=True)
            continue

        key = ObservedCombination.get_key(seq_regions, group=read_group)
        if not key in counts:
            counts[key] = ObservedCombination(seq_regions, group=read_group)
        counts[key].count += 1

    return list(counts.values()), filtered_seqs

def compare_to_library(counts, library, threshold=3, threads=1, progress=None, chunksize=None):
    """Determine the closest library sequences for each observed combination"""
    progress = _get_progress(progress)
    # Use multiprocessing heuristic for default chunksize
    if chunksize is None:
        chunksize, extra = divmod(len(counts), threads * 4)
        chunksize += bool(extra)

    if threads == 1 or len(counts) < chunksize:
        for i in progress(counts, desc="Library comparison"):
            i.compare_library(library, threshold=threshold)
    else:
        with multiprocessing.Pool(threads) as pool:
            # Use imap to allow tqdm progress bar

            f = _MultithreadCompareWrapper(library, threshold)
            count_iter = pool.imap(f, counts, chunksize=chunksize)
            counts = list(progress(count_iter, total=len(counts), desc="Library comparison"))

    return counts

def to_tsv_str(x, float_precision=5):
    """
    Convert x into a TSV compatible string
    """
    if isinstance(x, (str, int)):
        return x
    if isinstance(x, float):
        return f"{x:.{float_precision}f}"
    if x is None:
        return "NA"
    if isinstance(x, (list, tuple)):
        return ",".join(x)
    warnings.warn(f"Unrecognised type ({type(x)}), using raw string coercion")
    return str(x)

def write_counts(counts, outfile=sys.stdout, sort=False, progress=False):
    """
    Write counts to outfile
    """
    progress = _get_progress(progress)
    if sort:
        counts = sorted(counts, key=lambda x: x.count, reverse=True)

    print("group", end="\t", file=outfile)
    region_names = [r.name for r in counts[0].regions.values()]
    for name in region_names:
        print(name, f"{name}_in_library", f"{name}_nearest", f"{name}_distance",
              sep="\t", end="\t", file=outfile)
    print("combination_in_library", "nearest_combination_in_library", "count",
          sep="\t", end="\n", file=outfile)

    for comb in progress(counts, desc="Writing TSV"):
        print(to_tsv_str(comb.group), end="\t", file=outfile)
        for name in region_names:
            reg = comb.regions[name]
            print(*map(to_tsv_str, (reg.sequence, reg.in_library, reg.nearest, reg.distance)),
                  sep="\t", end="\t", file=outfile)
        print(*map(to_tsv_str, (comb.in_library, comb.nearest_in_library, comb.count)),
              sep="\t", end="\n", file=outfile)

def create_default_tqdm(newline=False, **kwargs):
    """Create TQDM factories with default args"""
    _tqdm = tqdm
    if newline:
        _tqdm.status_printer = staticmethod(_file_status_printer)

    outer_args = kwargs

    def f(*args, **kwargs):
        kwargs = {**outer_args, **kwargs}
        return _tqdm(*args, **kwargs)

    return f

def _get_progress(progress):
    """
    Get default TQDM progress bars.

    None = the standard bar
    False = disabled bar
    Otherwise passed through
    """
    if progress is None:
        progress = tqdm

    if not progress:
        progress = create_default_tqdm(disable=True)

    return progress

def _file_status_printer(file):
    """
    Patched TQDM status printer that uses \\n instead of \\r
    """
    fp = file
    fp_flush = getattr(fp, 'flush', lambda: None)  # pragma: no cover
    if fp in (sys.stderr, sys.stdout):
        getattr(sys.stderr, 'flush', lambda: None)()
        getattr(sys.stdout, 'flush', lambda: None)()

    def fp_write(s):
        fp.write(_unicode(s))
        fp_flush()

    last_len = [0]

    def print_status(s):
        len_s = disp_len(s)
        fp_write('\n' + s + (' ' * max(last_len[0] - len_s, 0)))
        last_len[0] = len_s

    return print_status

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