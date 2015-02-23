#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Uniquify the sequences. 
This only exists because 'uniq' doesn't seem to have a flexible field separator
"""

import os, logging
import subprocess


def uniquify(input_filename, output_filename):
    """
    Sort and uniquify, but uniquify only on the sequence part
    """

    command = ['sort', input_filename]
    process = subprocess.Popen(command, close_fds=True, stdout=subprocess.PIPE, bufsize=16384)

    last_sequence = None
    with open(output_filename,'w') as output_file:
        for line in process.stdout:
            sequence, description = line.split(',', 1)
            if sequence != last_sequence:
                last_sequence = sequence
                output_file.write(line)

    sort_exit_code = process.wait()
    if sort_exit_code != 0:
        logging.error("sort exited with %s" % sort_exit_code)
    return sort_exit_code

def main(args):
    import argparse
    import logging
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', dest="verbosity", action="count", default=0,
        help="How verbose logging should be. The more invocations, the more verbose (up to three)")
    parser.add_argument('-o', '--output', dest="output_filename", required=True,
        help="""Filename to save the output to""")    
    parser.add_argument("sequence", nargs=1,
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

    return uniquify(parameters.sequence[0], parameters.output_filename)


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))