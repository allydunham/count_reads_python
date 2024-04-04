# Python Structured Read Counter

Python based structured sequence read counter split out from Swallow into it's own package

## Installation

``` bash
git clone gitlab.internal.sanger.ac.uk:ad44/count-reads-python.git
cd count_reads_python
pip install .
```

## Usage

``` bash
count_reads --help

usage: count_reads [-h] [--regions REGIONS] [--library LIBRARY] [--group GROUP] [--sort] [--min_phred MIN_PHRED] [--max_distance MAX_DISTANCE] [--format {auto,fasta,fastq}]
                   [--gzip {auto,gzip,none}] [--threads THREADS] [--chunksize CHUNKSIZE] [--verbose] [--debug]
                   F1 [F2]

Generate count tables from sequence data. Read (paired) fasta/fastq files and count the occurance of each read or sub-region thereof. Gzip compression is supported and filetype is
detected via the file extension. Observed reads can be compared to an expected library, including fuzzy matching. Outputs a TSV count table. Regions are defined in a TSV with four
columns. Lines starting with # are skipped: Name - The region name Read - Match the reverse (entry is any leading substring of reverse) or forward (otherwise) read Start - 0-based
inclusive start of the region End - 0-based exclusive end of the read The library to compare to is also defined in a TSV file, with one column per region to be compared to the library.
Each row contains the sequences of these regions from one oligo (the combination observed in reads will also be checked). Optional per region max accepted distances can be given at the
end of the headers after a bar (|). Output can be grouped by features of the sequence name by supplying a RegEx with a capture group to the --group argument. The result of the capture
group is used to group output counts as if it were an additional region.

positional arguments:
  F1                    File to count sequences from
  F2                    File to count sequences from (default: None)

options:
  -h, --help            show this help message and exit
  --regions REGIONS, -r REGIONS
                        Path to TSV file defining regions of interest (default: None)
  --library LIBRARY, -l LIBRARY
                        Path to TSV file defining expected library (default: None)
  --group GROUP, -g GROUP
                        Group output by this feature (default: None)
  --sort, -s            Sort output by descending count (default: False)
  --min_phred MIN_PHRED, -p MIN_PHRED
                        Filter reads with mean Phred score below this threshold (default: 0)
  --max_distance MAX_DISTANCE, -d MAX_DISTANCE
                        Only assign nearest sequences with at most this many substitutions (default: 3)
  --format {auto,fasta,fastq}, -f {auto,fasta,fastq}
                        Input file(s) format (default: auto)
  --gzip {auto,gzip,none}, -z {auto,gzip,none}
                        Input file(s) is gzipped (default: auto)
  --threads THREADS, -t THREADS
                        Number of threads to use (default: 1)
  --chunksize CHUNKSIZE, -c CHUNKSIZE
                        Chunksize for parallel processing (default: None)
  --verbose, -v         Verbose logging (default: False)
  --debug               Debug mode (default: False)
```

## Todo List

- Add docs and better explanation of library/regions format
- Tests
- Optimise (probably convert to Rust)
- Add index column to libraries and include in output
- Allow multiple input files
- Options to collapse and summarise results (e.g. nearest hits, unique only, only in library). Possibly separate script
- Summary statistics (mapping quality, recombination etc.). Probably separate script
