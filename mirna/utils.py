"""
Miscellanea
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import logging
import re

def convertToRNA(sequence):
    return sequence.upper().replace('T','U')

def normalisedRNA(file):
    return ''.join(convertToRNA(line) for line, _ in extractSequences(file))


def dnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'T' or base1 == 'T' and base2 == 'A'

def rnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'U' or base1 == 'U' and base2 == 'A'


def generateRNACombinations(length):
    """
    Generate all sequences of ACGU length long
    """
    if length == 0:
        return []

    sequences = ['A', 'C', 'U', 'G']

    for _ in range(2, length+1):
        new_sequences = []
        for sequence in sequences:
            new_sequences.extend('%s%s' % (sequence, base) for base in 'ACUG')
        sequences = new_sequences

    return sequences


def flexibleOpen(filename):
    """
    Open a file, decompressing it if it looks to be compressed
    """

    if filename.endswith('.gz'):
        import gzip
        file = gzip.open(filename, 'r')
    elif filename.endswith('.bz2'):
        import bz2
        file = bz2.BZ2File(filename, 'r', 16384)
    else:
        file = open(filename)

    return file



def extractSequences(file):
    """
    Iterator over fasta or fastq files. Generates pairs of space-less lines of contig, descriptions (or None in straight data files)
    """
    fastq = all_data = None # we don't know whether it is a fastq or a fasta yet
    spacer = re.compile(r'\s')

    def despace(line):
        return spacer.sub('', line)

    # find out if it's fasta or fastq or something else

    first_lines = []

    # grab two lines with content
    try:
        while len(first_lines) < 2:
            line = next(file)
            if line.strip():
                first_lines.append(line)
    except StopIteration:
        # we didn't even get to two
        pass

    count = len(first_lines)
    if count == 0:
        # nothing there
        return
    elif count == 1:
        # one line. probably just a data line
        if not first_lines[0].startswith('>') and not first_lines[0].startswith('@'):
            contig = despace(first_lines[0])
            if contig and contig[0] in 'ACGTUacgtu':
                yield contig, None
        return
    else:
        if not first_lines[0].startswith('>') and not first_lines[0].startswith('@'):
            # more than one line, and first line is not a description. Going to check for all data
            # otherwise bail
            all_data = True
            for line in first_lines:
                contig = despace(line)
                if contig:
                    if contig[0] in 'ACGTUacgtu':
                        yield contig, None
                    else:
                        all_data = False
                        logging.error("Don't recognise format of file")
                        return
        else:
            all_data = False
            if first_lines[0].startswith('>'):
                # fasta
                fastq = False
                contigs = [despace(first_lines[1])]
                if contigs[0]:
                    # consume until the next line
                    try:
                        while True:
                            next_part = next(file)
                            if next_part.startswith('>'):
                                # end of first record, but keep the description
                                description = next_part[1:].strip()
                                yield ''.join(contigs), first_lines[0][1:].strip()
                                break
                            else:
                                contigs.append(despace(next_part))
                    except StopIteration:
                        # end of file. Only one record
                        yield ''.join(contigs), first_lines[0][1:].strip()
                        return
                else:
                    # we were all set for something standard, but no
                    logging.error("Don't recognise format of file")
                    return
            else:
                # fastq
                fastq = True
                contig = despace(first_lines[1])
                if contig:
                    yield contig, first_lines[0][1:].strip()
                else:
                    # we were all set for something standard, but no
                    logging.error("Don't recognise format of file")
                    return

                # consume the rest of the first record                    
                try:
                    separator = next(file)                
                    quality = next(file)
                except StopIteration:
                    # that's an early stop
                    return

    if all_data:
        for line in file:
            contig = despace(line)
            if contig:
                yield contig, None
    elif not fastq:
        # fasta
        contigs = []
        for line in file:
            if line.startswith('>'):
                if contigs:
                    yield ''.join(contigs), description
                description = line[1:].strip()
                contigs = []
            else:
                contigs.append(despace(line))
        if contigs:
            # leftovers
            yield ''.join(contigs), description
    else:
        # fastq
        try:
            while True:
                description = next(file)
                sequence = next(file)
                separator = next(file)
                quality = next(file)

                contig = despace(sequence)
                yield contig, description[1:].strip()
        except StopIteration:
            pass


