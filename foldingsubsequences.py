#!/usr/bin/env python
#coding: utf8

"""
Find substrings of each sequence that have ends that match (Watson-Crick match that is)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys, os
import logging
import re

from utils import rnaCrickMatch, extractSequences, flexibleOpen, convertToRNA, generateRNACombinations

DefaultEndBasePairs = 3
DefaultMinLength = 60
DefaultMaxLength = 114

NonBaseFinder = re.compile(r'[^ACUG]')

def findSubstringsThatMatch(line, end_base_pairs, min_length, max_length, max_matches=None):
    matches = []

    # Find the continuous sequences of bases, so that we don't have to check in the inner loop
    clean_parts = NonBaseFinder.split(line)

    for sequence in clean_parts:
        if not sequence:
            continue
        sequence_length = len(sequence)
        for start in range(sequence_length - min_length + 1):
            max_available_length = min(sequence_length - start, max_length)

            # go from largest to smallest to bias big in case we are only looking for a couple
            for length in range(max_available_length, min_length-1, -1):
                substring = line[start:start+length]

                start_match = substring[:end_base_pairs]
                end_match = substring[-end_base_pairs:][::-1]

                if all(rnaCrickMatch(base1, base2) for base1, base2 in zip(start_match, end_match)):
                    matches.append(substring)
                    if max_matches and len(matches) >= max_matches:
                        return matches
    return matches


def generatePossibleSubsequences(input_filename, output_filename, end_base_pairs, min_length, max_length, all_combinations=True, split_level=0):
    """
    Emit all substrings of between minLength and maxLength characters from lines in filenames 
    where the endBasePairs pairs at either end match with each other (Watson-Crick style)
    """

    parsed_sequences = set()

    try:
        fasta_file = flexibleOpen(input_filename)
    except (OSError, IOError) as error:
        logging.error("Could not open %s" % input_filename)
        return False

    if all_combinations:
        max_matches = None
    else:
        max_matches = 1

    for line in extractSequences(fasta_file):
        line = convertToRNA(line)

        for result in findSubstringsThatMatch(line, end_base_pairs, min_length, max_length, max_matches=max_matches):
            parsed_sequences.add(result)

    fasta_file.close()

    if not parsed_sequences:
        return True

    if split_level == 0:
        try:
            with open(output_filename,'w') as output_file:
                for sequence in parsed_sequences:
                    output_file.write("%s\n" % sequence)
        except (IOError, OSError) as error:
            logging.error("Problems trying to create %s: %s" % (output_filename, error))
            return False
    else:

        # Split the sequences based on their split_level starting bases

        file_mapper = {}
        prefices = generateRNACombinations(split_level)

        for prefix in prefices:
            output_file = open("%s.%s" % (output_filename, prefix), 'w')
            file_mapper[prefix] = output_file


        for sequence in parsed_sequences:
            file_mapper[sequence[:split_level]].write("%s\n" % sequence)


        for output_file in file_mapper.values():
            output_file.close()

    return True


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
