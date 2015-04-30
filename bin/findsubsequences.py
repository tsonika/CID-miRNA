#!/usr/bin/env python
#coding: utf8

"""
Find substrings of each sequence that have ends that match (Watson-Crick match that is)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys, os
import logging

from mirna.utils import rnaCrickMatch, extractSequences, flexibleOpen, convertToRNA, generateRNACombinations
from mirna.foldingsubsequences import generatePossibleSubsequences, DefaultMinLength, DefaultMaxLength, DefaultEndBasePairs


def main(args):

    from argparse import ArgumentParser
    
    def levelFromVerbosity(verbosity):
        if verbosity == 0:
            level = logging.ERROR
        elif verbosity == 1:
            level = logging.WARNING
        elif verbosity == 2:
            level = logging.INFO
        else:
            level = logging.DEBUG

        return level
    
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", dest="verbosity", default=0, action="count",
                      help="Verbosity.  Invoke many times for higher verbosity")
    parser.add_argument('-o', '--output', dest="output_filename", required=True,
        help="""Filename to save the output to""")    
    parser.add_argument("-b", "--end-base-pairs", dest="end_base_pairs", type=int, default=DefaultEndBasePairs,
        help="Number of bases to force pair at the ends (default: %(default)s)")    
    parser.add_argument("-1", "--one-candidate", dest="one_candidate_per_line", default=False,
        action="store_true",
        help="Generate only one candidate miRNA per sequence line instead of all possible subsequences")    
    parser.add_argument("-m", "--min-length", dest="min_length", type=int, default=DefaultMinLength,
        help="Minimum length allowed (default: %(default)s)")
    parser.add_argument("-x", "--max-length", dest="max_length", type=int, default=DefaultMaxLength,
        help="Maximum length allowed (default: %(default)s)")
    parser.add_argument("-s", "--split", dest="split_level", type=int, default=0,
        help="How many levels to split the output files into (default: %(default)s)")
    parser.add_argument("sequence", nargs=1,
                      help="Filename containing FASTA or FASTQ sequences")


    parameters = parser.parse_args(args)
    # logger doubles as a section timer
    logging.basicConfig(level=levelFromVerbosity(parameters.verbosity), format="%(asctime)s:%(levelname)s:%(name)s:%(message)s")

    if generatePossibleSubsequences(parameters.sequence[0], parameters.output_filename, 
        end_base_pairs=parameters.end_base_pairs, min_length=parameters.min_length, 
        max_length=parameters.max_length, all_combinations=not parameters.one_candidate_per_line, 
        split_level=parameters.split_level):
        return 0
    else:
        return 1

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
