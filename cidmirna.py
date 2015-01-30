#!/usr/bin/env python
#coding: utf8

"""
Extract candidate miRNA sequences from fasta or fastq sequences
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys, os
import logging
import tempfile
import shutil
import subprocess
import re
import time


DefaultWindowLength = 125
DefaultUpperDGCutoff = -13.559
DefaultWindowStepSize = 10
DefaultEndBasePairs = 3
DefaultMinLength = 60
DefaultMaxLength = 114
DefaultGrammarScoreCutoff = -0.609999
DefaultStructuralScoreCutoff = 23
DefaultProbabilitiesFilename = "CFGprobabilities.txt"

# We are packaging old binaries until we recreate the source
NeedPreload = ['newcyk2','vienna2struct']

#--------------------------------------------------------------------------------------------------------------------------
# This program calls the following six programs in succession for a full genome scan
#--------------------------------------------------------------------------------------------------------------------------
# Program Name          Input Arguments                         Needed in this program
#--------------------------------------------------------------------------------------------------------------------------
#[1] autoparseunique.pl
#               - Input file with only sequences to be parsed               N
#               - Number of bases to force pair at the ends             default = 3
#               - Lowest length allowed                         default = 60
#               - Highest length allowed                        default = 114
#
# [2] newcyk2
#               - Input file with only sequences to be run through the grammar      N
#               - Input probability file name                       N
#               - Output result filename                        N
#
# [3] cutoffpassscore
#               - Input filename with sequences and scores              N
#               - Cutoff grammar score                          default = -0.609999
#               - Output filename                           N
#
# [4] findlocus.pl
#               - Input filename from CutOff selector                   N
#               - Source fasta format filename                      N
#               - Output filename                           N
# [5] orderpos.pl
#               - Input filename with unordered positions               N
#
# [6] remove_overlap
#               - Input filename with overlaps                      N
#               - Output filename without overlaps                  N
# [7] convert2fasta.pl
#               - Input filename in the final results format                N
#
# [8] grammar2rfold.pl
#               - Input filename in the final results format                N
#
# [9] RNAfold
#               - Some parameters                           N (supplied here)
#               - Input filename with sequences and structural constraints      N
#
# [10] gawk
#               - Some parameters                           N
#               - Input filename with the Rfold structures              N
#
# [11] mergeLoops.pl
#               - Input filename with the Rfold structures              N
#
# [12] vienna2struct
#               - Input filename with the Rfold structures (merged loops)       N
#               - Output filename for storing the drawn structures          N
#
# [13] Scores4mStruct
#               - Input filename with the drawn structures              N
#               - Output filename for storing the scores                N
#
# [14] StructScorePass
#               - Input filename with the scores and structures             N
#               - Output filename for storing the final results             N
#
#--------------------------------------------------------------------------------------------------------------------------


def extractFasta(file):
    """
    Iterator over fasta files. Generates space-less lines of contig
    """
    spacer = re.compile(r'\s')
    for line in file:
        if not line.startswith('>'):
            contig = spacer.sub('', line)
            if contig:
                yield contig

def convertToRNA(sequence):
    return sequence.upper().replace('T','U')

def normalisedRNA(file):
    return ''.join(convertToRNA(line) for line in extractFasta(file))


def runCommand(command, filename, parameters, output_extension, output_is_stdout=False, input=None, manual_parameters=False, local=True):

    ld_library_path = None
    if local:
        script_path = os.path.dirname(__file__)
        if not script_path:
            script_path = '.'

        if command in NeedPreload:
            ld_library_path = script_path
        command = os.path.join(script_path, command)


    output_filename = "%s.%s" % (filename, output_extension)
    if manual_parameters:
        full_command = [command] + parameters
    else:
        full_command = [command, filename] + parameters

    if output_is_stdout:
        output_file = open(output_filename, 'w')
        output = output_file
    else:
        full_command.append(output_filename)
        output = None

    environment = dict(os.environ)
    if ld_library_path:
        environment['LD_LIBRARY_PATH'] = "%s%s" % (ld_library_path, ':%s' % environment.get['LD_LIBRARY_PATH'] if "LD_LIBRARY_PATH" in environment else '')
    try:
        logging.info("Running %s" % " ".join(full_command))
        subprocess.check_call(full_command, stdout=output, stdin=input, close_fds=True, env=environment)
        logging.info("Finished running %s" % " ".join(full_command))
    except subprocess.CalledProcessError as error:
        logging.error("%s exited with code %s" % (' '.join(full_command), error.returncode))
        return False

    if output_is_stdout:
        output_file.close()

    return output_filename



def rnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'U' or base1 == 'U' and base2 == 'A'

def dnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'T' or base1 == 'T' and base2 == 'A'


def findSubstringsThatMatch(line, endBasePairs, minLength, maxLength):
    matches = []
    line_length = len(line)
    for start in range(line_length - minLength + 1):
        maxAvailableLength = min(line_length - start, maxLength)
        for length in range(minLength, maxAvailableLength + 1):
            substring = line[start:start+length]
            if all(character in 'ACUG' for character in substring):
                startMatch = substring[:endBasePairs]
                endMatch = substring[-endBasePairs:][::-1]

                if all(rnaCrickMatch(base1, base2) for base1, base2 in zip(startMatch, endMatch)):
                    matches.append(substring)
            else:
                # if all characters are not legal now, they won't be when we extend the string
                break


    return matches



def parser4auto(filename, endBasePairs, minLength, maxLength):
    """
    Python version of parser4auto.cpp. 
    Emit all substrings of between minLength and maxLength characters from lines in filename 
    where the endBasePairs pairs at either end match with each other (Watson-Crick style)
    """

    MaxListSize = 850

    parsed_sequences = set()
 
    output_counter = 1

    try:
        fasta_file = open(filename)
    except (OSError, IOError) as error:
        logging.error("Could not open %s" % filename)
        return False

    for line in extractFasta(fasta_file):
        logging.debug("line: %s" % line)

        line = convertToRNA(line)
        for result in findSubstringsThatMatch(line, endBasePairs, minLength, maxLength):
            parsed_sequences.add(result)
                
            if len(parsed_sequences) > MaxListSize:
                break

        if len(parsed_sequences) > MaxListSize:
            break


    fasta_file.close()

    if not parsed_sequences:
        return False

    output_filename = "%s.%s" % (filename, output_counter)
    try:
        with open(output_filename,'w') as output_file:
            for sequence in parsed_sequences:
                output_file.write("%s\n" % sequence)
    except (IOError, OSError) as error:
        logging.error("Problems trying to create %s: %s" % (output_filename, error))
        return False

    return output_filename



def autoParseUnique(filename, endBasePairs, minLength, maxLength):
    """
    Automate the process of parsing sequences from the parser
    """


    MaxListSize = 850

    parsed_sequences = set()
 
    output_counter = 1

    try:
        fasta_file = open(filename)
    except (OSError, IOError) as error:
        logging.error("Could not open %s" % filename)
        return False

    current_directory = os.path.abspath(os.curdir)

    temp_directory = tempfile.mkdtemp(prefix='cidmirna')
    logging.debug("Temp directory: %s" % temp_directory)
    ParsedMarker = 'PARSED SEQUENCE:'
    for line in extractFasta(fasta_file):
        logging.debug("line: %s" % line)
        search_filename = os.path.join(temp_directory, "TempFile")
        search_file = open(search_filename, 'w')
        search_file.write("%s\n" % line)
        search_file.close()

        command = [os.path.join(current_directory, 'parser4auto'), search_filename, str(endBasePairs), str(minLength), str(maxLength)]


        try:
            logging.info("Running %s" % " ".join(command))
            output = subprocess.check_output(command, cwd=temp_directory, close_fds=True)
            for output_line in output.splitlines():
                if output_line.startswith(ParsedMarker):
                    rest = output_line[len(ParsedMarker):]
                    parsed_sequences.add(rest)
                    if len(parsed_sequences) > MaxListSize:
                        break

        except subprocess.CalledProcessError as error:
            logging.error("%s exited with code %s" % (' '.join(command), error.returncode))
            fasta_file.close()
            return False

        else:
            if len(parsed_sequences) > MaxListSize:
                break

    fasta_file.close()
    shutil.rmtree(temp_directory)


    output_filename = "%s.%s" % (filename, output_counter)
    try:
        with open(output_filename,'w') as output_file:
            for sequence in parsed_sequences:
                output_file.write("%s\n" % sequence)
    except (IOError, OSError) as error:
        logging.error("Problems trying to create %s: %s" % (output_filename, error))
        return False

    return output_filename


def runNewcyk(filename, probabilities_filename):
    return runCommand('newcyk2', filename, [probabilities_filename], 'grm')

def runCutoffPassscore(filename, score):
    return runCommand('cutoffpassscore', filename, [str(score)], "cof")


def findLoci(filename, fasta_filename):
    """
    Look for the given sequences returned from the CutOff selector program in the source Contig file, and
    writes another file that can serve as input for the Overlap Remover
    """

    output_filename = "%s.%s" % (filename, "loci")

    input_file = open(filename)

    contig_file = open(fasta_filename)
    rna = normalisedRNA(contig_file)
    contig_file.close()

    sequences = []

    #Run through the Input File and look for its position in the ContigFile
    for line in input_file:
        if line.startswith('Sequence :'):
            sequence_line = next(input_file).strip()
            if sequence_line:
                score_line = next(input_file)
                score_line = score_line.strip()

                _, _, _, sequence_length, _, _, _, normalised_sequence_score, _, _, sequence_score = score_line.split()
                sequence_line = sequence_line.upper()

                # Find all the locations of the sequence in the original file
                found = False
                last_index = 0
                while last_index >= 0:
                    last_index = rna.find(sequence_line, last_index)
                    if last_index >= 0:
                        found = True

                        sequences.append((last_index, """ Predicted miRNA\tPosition: %s   Length = %s\n Score = %s\tNormalized Score = %s\n%s\n""" % (last_index, sequence_length, sequence_score, normalised_sequence_score, convertToRNA(sequence_line))))
                        last_index += 1

                if not found:
                    logging.error("A sequence in the input file was not found in the Source fasta.\n%s\nContinuing with other sequences...\n\n" % sequence_line)

    input_file.close()

    # sort the sequences by location
    sequences.sort()

    output_file = open(output_filename, "w")
    for sequence in sequences:
        output_file.write(sequence[1])
    output_file.close()    

    return output_filename


def convertToFasta(filename):
    output_filename = "%s.fasta" % filename

    input_file = open(filename)
    output_file = open(output_filename, 'w')

    for line in input_file:
        if 'Sequence' in line:
            sequence = next(input_file)

            score_line = next(input_file)
            bits = score_line.strip().split()
            length = bits[2]
            score = bits[9]
            normalised_score = bits[6]

            fasta_line = ">Predicted\tLength: %s\tScore: %s\tNormalised Score: %s\n%s" % (length, score, normalised_score, sequence)
            output_file.write(fasta_line)
    output_file.close()
    input_file.close()

    return output_filename


def grammarToRFold(filename):
    """
    Convert FASTA to something that is recognised by RNAfold with constraints options

    FASTA should be in the following format:
    >Predicted      Length: 100     Score: -60.1264 Normalised Score: -0.601264
    CCTTCTAGTGGCAAGAGTGACGTAAGTGATATGCGGAAATTTCTTTCCAAGCCTGCTTGGAGAAGCTTCCTCTGCCTGCTTCTCTTTGGCCACCTCCAGG
    >Predicted      Length: 75      Score: -45.5868 Normalised Score: -0.607824
    GTTTTCCCTCTTATGTCCAGCAAATGCTGCATGGAGCCCTGGAATTCTATGTGGAAAGCTAGGAAGAGGGAGAGC
    >Predicted      Length: 62      Score: -37.5368 Normalised Score: -0.605432
    GGGTCTTTGTGTCAATCTGAGCTCTGATGTCCACCTAGAGATTGGGTATCCACCTAAGGCCC

    """

    output_filename = "%s.4rfold" % filename

    input_file = open(filename)
    output_file = open(output_filename,'w')

    for line in input_file:
        stripped_line = line.strip()
        if not stripped_line or stripped_line.startswith("#"):
            continue

        if stripped_line.startswith('>'):
            bits = stripped_line.split()
            score = bits[4]
            normalised_score = bits[-1]

            output_file.write('>Score(%s)(%s)\n' % (score, normalised_score)) 
        else:
            output_file.write(line)
            line_length = len(line) - 1
            half_length = int(line_length / 2)
            output_file.write("%s%s\n" % ('<' * half_length, '.' * (line_length - half_length)))

    output_file.close()
    input_file.close()

    return output_filename


def runRNAFold(filename):
    input = open(filename)
    output_filename = runCommand('RNAfold', filename, ['-C', '--noPS'], 'rfold', output_is_stdout=True, manual_parameters=True, local=False, input=input)
    input.close()

    return output_filename


def removeMeanFreeEnergyValues(filename):
    """
    Remove the scores produced by RNAfold
    """
    return runCommand('gawk', filename, ['{print $1}', filename], 'nodg', output_is_stdout=True, manual_parameters=True, local=False)



def keepOneLoop(structure):
    """
    Keep only the centre loop in a structure. 
    """

    loop_centre = re.compile(r'\(\.*\)')

    
    loop_centres = []
    halfway = len(structure) / 2
    for match in loop_centre.finditer(structure):
        loop_centre = (match.start() + match.end()) / 2
        loop_centres.append((abs(halfway - loop_centre), loop_centre))

    closest_to_centre = min(loop_centres)
    centre_loop = closest_to_centre[1]

    # We've decided on the winner, now get rid of every other loop

    def straightenLoopsBefore(structure, end, closes_loop, opens_loop):
        fixed_until = end - 1
        while fixed_until >= 0:
            closing_loop = structure.rfind(closes_loop, 0, fixed_until)
            if closing_loop < 0:
                # done
                break

            closed = 1
            index = closing_loop - 1

            while index >= 0:
                if structure[index] == closes_loop:
                    closed += 1
                elif structure[index] == opens_loop:
                    closed -= 1
                    if closed == 0:
                        # all found, replace everything in between with dots
                        structure = structure[:index] + ('.' * (closing_loop - index + 1)) + structure[closing_loop+1:]
                        fixed_until = index - 1
                        break
                index -= 1
        return structure

    # First backwards
    structure = straightenLoopsBefore(structure, centre_loop, ')', '(')
    # Now forwards, by reversing the structure
    structure = structure[::-1]
    structure = straightenLoopsBefore(structure, centre_loop, '(', ')')
    structure = structure[::-1]

    return structure


def mergeLoops(filename):
    """
    Make sure there's only one loop per structure. If there are more, keep the centre one and 
    'straighten' the other ones out
    """

    output_filename = "%s.mloops" % filename

    divergent_loops = re.compile(r'\).*.\(')

    output_file = open(output_filename,'w')
    input_file = open(filename)
    for line in input_file:
        stripped_line = line.strip()
        if not line or stripped_line.startswith('#'):
            continue

        if stripped_line.startswith('>'):
            output_file.write(line)
        elif stripped_line[0] not in 'ACGTUacgtu':
            count = 0

            for _ in divergent_loops.finditer(stripped_line):
                count += 1
                break

            if count == 0:
                # only one loop. Write it as is
                output_file.write(line)
            else:
                # Multiple loops. There must be only one
                straightened_structure = keepOneLoop(stripped_line)
                output_file.write("%s\n" % straightened_structure)

        else:
            output_file.write(line)


    output_file.close()
    input_file.close()

    return output_filename



def filterOnScores(filename, cutoff_score):
    """
    Filter all structures below the cutoff_score out
    """

    output_filename = "%s.pass" % filename

    input_file = open(filename)
    output_file = open(output_filename, 'w')

    for line in input_file:
        if line.startswith('>'):
            sequence = next(input_file)
            struct_lines = []
            for _ in range(5):
                struct_lines.append(next(input_file))

            blank_line = next(input_file).strip()
            if blank_line:
                struct_score_line = blank_line
            else:
                struct_score_line = next(input_file).strip()

            score = 0
            bits = struct_score_line.split()
            if bits[0] == 'SCORE:':
                score = int(bits[1])

            if score >= cutoff_score:
                output_file.write(line)
                output_file.write(sequence)
                for struct_line in struct_lines:
                    output_file.write(struct_line)
                output_file.write('%s\n' % struct_score_line)
    output_file.close()
    input_file.close()

    return output_filename


def structuresToFasta(filename):
    output_filename = "%s.fasta" % filename

    input_file = open(filename)
    output_file = open(output_filename, 'w')

    for line in input_file:
        if line.startswith('>'):
            parts = re.split(r'[>()\s+]', line)
            locus = parts[1]
            score = parts[2]
            normalised_score = parts[-1]

            sequence = next(input_file).replace('U','T')
            structure_lines = []
            for _ in range(5):
                structure_lines.append(next(input_file))

            structure_score_line = next(input_file)
            bits = structure_score_line.split()
            if bits[0] == 'SCORE:':
                structure_score = bits[1].strip()
            else:
                structure_score = '0'
            
            fasta_line = ">Predicted-miRNA Position:%(locus)s Length:%(length)s GrammarScore:%(score)s NormalisedGrammarScore:%(normalise_score)s StructuralScore:%(structure_score)s\n" % {
                'locus' : locus,
                'length' : len(sequence.strip()),
                'score' : score,
                'normalise_score' : normalised_score,
                'structure_score' : structure_score
            }

            output_file.write(fasta_line)
            output_file.write(sequence)

            for structure_line in structure_lines:
                output_file.write(structure_line)

    output_file.close()
    input_file.close()

    return output_filename




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
    parser.add_argument("-t", "--train", dest="train", default=False, action="store_true",
                      help="Train model")
    parser.add_argument("-l", "--window-length", dest="windowLength", default=DefaultWindowLength, type=int,
        help="Fixed window length to scan with (default: %(default)s)")
    parser.add_argument("--window-step-size", dest="windowStepSize", type=int, default=DefaultWindowStepSize,
        help="Step size to use for shifting window (default: %(default)s)")
    parser.add_argument("-dG", "--upper-dg", dest="upperDGCutoff", type=float, default=DefaultUpperDGCutoff,
        help="Upper dG cutoff value (default: %(default)s)")
    parser.add_argument("-b", "--end-base-pairs", dest="endBasePairs", type=int, default=DefaultEndBasePairs,
        help="Number of bases to force pair at the ends (default: %(default)s)")
    parser.add_argument("-m", "--min-length", dest="minLength", type=int, default=DefaultMinLength,
        help="Minimum length allowed (default: %(default)s)")
    parser.add_argument("-x", "--max-length", dest="maxLength", type=int, default=DefaultMaxLength,
        help="Maximum length allowed (default: %(default)s)")
    parser.add_argument("-g", "--grammar-score", dest="grammarCutoff", type=float, default=DefaultGrammarScoreCutoff,
        help="Cutoff grammar score (default: %(default)s)")
    parser.add_argument("-s", "--structural-score", dest="structuralCutoff", type=float, 
        default=DefaultStructuralScoreCutoff,
        help="Cutoff structural score (default: %(default)s)")    
    parser.add_argument("-e", "--email", dest="email", default=None,
        help="Email address to send results to")    
    #parser.add_argument("--loci-filename", dest="lociFilename", default=None, required=True,
    #    help="Name of file within which the loci need to be located")
    parser.add_argument("--probabilities-filename", dest="probabilitiesFilename", default=DefaultProbabilitiesFilename,
        help="Filename where to find the CFG probabilities (default: %(default)s)")
    parser.add_argument("sequences", nargs=1,
                      help="Filenames containing FASTA or FASTQ sequences")


    parameters = parser.parse_args(args)
    # logger doubles as a section timer
    logging.basicConfig(level=levelFromVerbosity(parameters.verbosity), format="%(asctime)s:%(levelname)s:%(name)s:%(message)s")

    result = 0

    filename = parameters.sequences[0]
    full_filename = os.path.expanduser(filename)
    if not os.path.isfile(full_filename):
        parser.error("%s is not a file" % filename)

    start_time = time.time()

    parsed_filename = parser4auto(full_filename, parameters.endBasePairs, parameters.minLength, parameters.maxLength)

    if not parsed_filename:
        logging.error("Couldn't parse sequences from %s." % filename)
        return 1

    grammar_filename = runNewcyk(parsed_filename, parameters.probabilitiesFilename)
    if not grammar_filename:
        logging.error("Grammar run failed on %s" % parsed_filename)
        return 1

    cutoff_filename = runCutoffPassscore(grammar_filename, parameters.grammarCutoff)
    if not cutoff_filename:
        logging.error("Cuttoff run failed on %s" % grammar_filename)
        return 1

    # loci_filename = findLoci(cutoff_filename, full_filename)
    # if not loci_filename:
    #     logging.error("Loci finder failed on %s" % cutoff_filename)
    #     return 1

    # no_overlap_filename = runCommand('remove_overlap', loci_filename, [], 'nooverlap')
    # if not no_overlap_filename:
    #     logging.error("Problems removing overlaps on %s" % loci_filename)
    #     return 1

    grammar_fasta_filename = convertToFasta(cutoff_filename)
    if not grammar_fasta_filename:
        logging.error("Problems converting %s to fasta" % cutoff_filename)
        return 1

    rfold_filename = grammarToRFold(grammar_fasta_filename)
    if not rfold_filename:
        logging.error("Problems converting FASTA to something for rfold on %s" % grammar_fasta_filename)
        return 1

    rfold_output_filename = runRNAFold(rfold_filename)
    if not rfold_output_filename:
        logging.error("Problems running RNAfold on %s" % rfold_filename)
        return 1

    nodg_filename = removeMeanFreeEnergyValues(rfold_output_filename)
    if not nodg_filename:
        logging.error("Problems removing mean free energy values on %s" % rfold_output_filename)
        return 1

    merged_loops_filename = mergeLoops(nodg_filename)
    if not merged_loops_filename:
        logging.error("Problems merging loops on %s" % nodg_filename)
        return 1

    vienna_filename = runCommand('vienna2struct', merged_loops_filename, [], 'struct')
    if not vienna_filename:
        logging.error("Problems running vienna2struct on %s" % merged_loops_filename)
        return 1


    scores_filename = runCommand('Scores4mStruct', vienna_filename, [], 'scores')
    if not scores_filename:
        logging.error("Problems running Scores4mStruct on %s" % vienna_filename)
        return 1


    filtered_filename = filterOnScores(scores_filename, parameters.structuralCutoff)
    if not filtered_filename:
        logging.error("Problems filtering structures on %s" % scores_filename)
        return 1

    fasta_filename = structuresToFasta(filtered_filename)
    if not fasta_filename:
        logging.error("Problems converting the structures to fasta files on %s" % filtered_filename)
        return 1


    final_fasta_filename = "%s.final.fasta" % full_filename
    os.rename(fasta_filename, final_fasta_filename)

    end_time = time.time()

    logging.info("Took %s seconds" % (end_time - start_time))

    print("FINAL OUTPUT FILES:\n\tPost-grammar cof: %s\n\tPost-structure cof: %s\n" % (grammar_fasta_filename, final_fasta_filename))


    return result

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
