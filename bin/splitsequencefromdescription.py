#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Split sequence and descriptions into separate files
"""


def split_sequences(input_filename, sequence_filename, description_filename):

    with open(input_filename) as input_file, open(sequence_filename, 'w') as sequence_file, \
    open(description_filename, 'w') as description_file:
        for line in input_file:
            sequence, description = line.split(',', 1)
            sequence_file.write("%s\n" % sequence)
            description_file.write(description)

    return 0

def main(args):
    import argparse
    import logging
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', dest="verbosity", action="count", default=0,
        help="How verbose logging should be. The more invocations, the more verbose (up to three)")
    parser.add_argument('-s', '--sequence', dest="sequence_filename", required=True,
        help="""Filename to save the sequences to""")    
    parser.add_argument('-d', '--descriptions', dest="description_filename", required=True,
        help="""Filename to save the descriptions to""")        
    parser.add_argument("input", nargs=1,
                      help="Filename containing sequences,description pairs")        
    parameters = parser.parse_args(args)

    def verbosityToLoggingLevel(verbosity):
        if verbosity <= 0:
            return logging.ERROR
        elif verbosity <= 1:
            return logging.WARNING
        elif verbosity <= 2:
            return logging.INFO
        else:
            return logging.DEBUG

    logging.basicConfig(level=verbosityToLoggingLevel(parameters.verbosity))


    return split_sequences(parameters.input[0], parameters.sequence_filename, parameters.description_filename)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))