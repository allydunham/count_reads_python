"""
A python CLI tool and library for counting and processing structured sequence reads
"""

# import logging

from count_reads.lib import *

__all__ = ["ReadPair", "SequenceParser", "SequenceLibrary", "Region", "ObservedRegion",
           "ObservedCombination", "FilterCount", "count_seqs", "compare_to_library",
           "write_counts"]

