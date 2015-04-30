#!/usr/bin/env python
#coding: utf8

"""
Find substrings of each sequence that have ends that match (Watson-Crick match that is)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys, os
import logging
import re

from mirna.utils import rnaCrickMatch, extractSequences, flexibleOpen, convertToRNA, generateRNACombinations

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

    try:
        fasta_file = flexibleOpen(input_filename)
    except (OSError, IOError) as error:
        logging.error("Could not open %s" % input_filename)
        return False

    if all_combinations:
        max_matches = None
    else:
        max_matches = 1


    # Open output files
    file_mapper = {}
    if split_level == 0:
        try:
            file_mapper[''] = open(output_filename,'w')
        except (IOError, OSError) as error:
            logging.error("Problems trying to create %s: %s" % (output_filename, error))
            return False
    else:

        # Split the sequences based on their split_level starting bases
        prefices = generateRNACombinations(split_level)

        for prefix in prefices:
            output_file = open("%s.%s" % (output_filename, prefix), 'w')
            file_mapper[prefix] = output_file


    for line, description in extractSequences(fasta_file):
        line = convertToRNA(line)
        if description is None:
            description = ''
        for result in findSubstringsThatMatch(line, end_base_pairs, min_length, max_length, max_matches=max_matches):
            file_mapper[result[:split_level]].write("%s,%s\n" % (result, description))

    fasta_file.close()

    for output_file in file_mapper.values():
        output_file.close()



    return True


