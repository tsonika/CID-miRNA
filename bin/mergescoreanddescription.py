#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Merge sequence and descriptions files into one file
"""


def merge(score_filename, description_filename, output_filename):

    with open(score_filename) as score_file, open(description_filename) as description_file, \
    open(output_filename, 'w') as output_file:
        for score_line in score_file:
            if score_line.strip().startswith("Sequence :"):
                try:
                    description = next(description_file)
                except StopIteration:
                    # newcyk2 sometimes repeats the last sequence. Harmless, just stop now
                    break
                sequence = next(score_file)
                score = next(score_file)
                output_file.write(" Sequence : %s" % description)
                output_file.write(sequence)
                output_file.write(score)

    return 0

def main(args):
    import argparse
    import logging
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', dest="verbosity", action="count", default=0,
        help="How verbose logging should be. The more invocations, the more verbose (up to three)")
    parser.add_argument('-s', '--score', dest="score_filename", required=True,
        help="""Filename where scores reside""")    
    parser.add_argument('-d', '--descriptions', dest="description_filename", required=True,
        help="""Filename where descriptions reside""")        
    parser.add_argument('-o', '--output', dest="output_filename", required=True,
        help="""Output filename""")            
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


    return merge(parameters.score_filename, parameters.description_filename, parameters.output_filename)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))