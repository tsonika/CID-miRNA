#!/usr/bin/env python
#coding: utf8

"""
Extract candidate miRNA sequences from fasta or fastq sequences
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys, os
import logging
import subprocess
import re
import time

from mirnastructure import diagramFromStructure

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
NeedPreload = ['newcyk2']

#--------------------------------------------------------------------------------------------------------------------------
# This program calls the following six programs in succession for a full genome scan
#--------------------------------------------------------------------------------------------------------------------------
# Program Name          Input Arguments                         Needed in this program
#--------------------------------------------------------------------------------------------------------------------------
# [1] newcyk2
#               - Input file with only sequences to be run through the grammar      N
#               - Input probability file name                       N
#               - Output result filename                        N
#
# [2] cutoffpassscore
#               - Input filename with sequences and scores              N
#               - Cutoff grammar score                          default = -0.609999
#               - Output filename                           N
#
#
# [3] RNAfold
#               - Some parameters                           N (supplied here)
#               - Input filename with sequences and structural constraints      N
#
# [4] gawk
#               - Some parameters                           N
#               - Input filename with the Rfold structures              N
#
#
# [5] Scores4mStruct
#               - Input filename with the drawn structures              N
#               - Output filename for storing the scores                N
#
# [6] StructScorePass
#               - Input filename with the scores and structures             N
#               - Output filename for storing the final results             N
#
#--------------------------------------------------------------------------------------------------------------------------

class StandardRunner(object):
    def run(self, command, output_filename=None, error_filename=None, input_filename=None, environment=None):
        parameters = {
            'close_fds' : True
        }

        if output_filename:
            output_file = open(output_filename, 'w')
            parameters['stdout'] = output_file
        if input_filename:
            input_file = open(input_filename)
            parameters['stdin'] = input_file
        if error_filename:
            if error_filename != output_filename:
                error_file = open(error_filename, 'w')
            else:
                error_file = output_file
            parameters['stderr'] = error_file

        if environment:
            parameters['env'] = environment

        logging.debug("Running %s" % " ".join(command))
        exit_code = subprocess.call(command, **parameters)
        logging.debug("Finished running %s with exit code %s" % (" ".join(command), exit_code))


        if input_filename:
            input_file.close()            
        if output_filename:
            output_file.close()
        if error_filename and error_filename != output_filename:
            error_file.close()

        return exit_code


    def max_processes(self):
        import multiprocessing
        return multiprocessing.cpu_count()



class Runner(object):

    Runner = StandardRunner()
    MaxProcesses = None
    OutputDirectory = '.'

    @classmethod
    def runCommand(cls, command, filename, parameters, output_extension, output_is_stdout=False, input_filename=None, manual_parameters=False, local=True):

        """
        Standard command run wrapper.
        """

        ld_library_path = None
        if local:
            script_path = os.path.dirname(__file__)
            if not script_path:
                script_path = '.'

            script_path = os.path.join(script_path, 'bin')

            if command in NeedPreload:
                ld_library_path = script_path
            command = os.path.join(script_path, command)



        output_filename = "%s.%s" % (filename, output_extension)
        if manual_parameters:
            full_command = [command] + parameters
        else:
            full_command = [command, filename] + parameters

        if output_is_stdout:
            piped_output = output_filename
        else:
            full_command.append(output_filename)
            piped_output = None

        environment = dict(os.environ)
        if ld_library_path:
            environment['LD_LIBRARY_PATH'] = "%s%s" % (ld_library_path, ':%s' % environment.get('LD_LIBRARY_PATH') if "LD_LIBRARY_PATH" in environment else '')

        exit_code = cls.Runner.run(full_command, output_filename=piped_output, 
            input_filename=input_filename, environment=environment)
        if exit_code != 0:
            return False

        return output_filename

    @classmethod
    def max_processes(cls):
        if cls.MaxProcesses:
            return cls.MaxProcesses
        else:
            return cls.Runner.max_processes()


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
    Iterator over fasta or fastq files. Generates space-less lines of contig
    """
    fastq = all_data = None # we don't know whether it is a fastq or a fasta yet
    spacer = re.compile(r'\s')

    def despace(line):
        return spacer.sub('', line)

    # find out if it's fasta or fastq or something else

    first_lines = []
    try:
        first_lines.append(next(file))
        first_lines.append(next(file))
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
                yield contig
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
                        yield contig
                    else:
                        all_data = False
                        logging.error("Don't recognise format of file")
                        return
        else:
            all_data = False
            if first_lines[0].startswith('>'):
                # fasta
                fastq = False
                contig = despace(first_lines[1])
                if contig:
                    yield contig
                else:
                    # we were all set for something standard, but no
                    logging.error("Don't recognise format of file")
                    return
            else:
                # fastq
                fastq = True
                contig = despace(first_lines[1])
                if contig:
                    yield contig                
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
                yield contig
    elif not fastq:
        # fasta
        for line in file:
            if not line.startswith('>'):
                contig = despace(line)
                if contig:
                    yield contig
    else:
        # fastq
        try:
            while True:
                description = next(file)
                sequence = next(file)
                separator = next(file)
                quality = next(file)

                contig = despace(sequence)
                yield contig
        except StopIteration:
            pass


def convertToRNA(sequence):
    return sequence.upper().replace('T','U')

def normalisedRNA(file):
    return ''.join(convertToRNA(line) for line in extractSequences(file))



def rnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'U' or base1 == 'U' and base2 == 'A'

def dnaCrickMatch(base1, base2):
    return base1 == 'C' and base2 == 'G' or base1 == 'G' and base2 == 'C' or base1 == 'A' and base2 == 'T' or base1 == 'T' and base2 == 'A'


def findSubstringsThatMatch(line, end_base_pairs, min_length, max_length, max_matches=None):
    matches = []
    line_length = len(line)
    for start in range(line_length - min_length + 1):
        max_available_length = min(line_length - start, max_length)
        for length in range(min_length, max_available_length + 1):
            substring = line[start:start+length]
            if all(character in 'ACUG' for character in substring):
                start_match = substring[:end_base_pairs]
                end_match = substring[-end_base_pairs:][::-1]

                if all(rnaCrickMatch(base1, base2) for base1, base2 in zip(start_match, end_match)):
                    matches.append(substring)
                    if max_matches and len(matches) >= max_matches:
                        return matches
            else:
                # if all characters are not legal now, they won't be when we extend the string
                break


    return matches



def generatePossibleSubsequences(filenames, end_base_pairs, min_length, max_length, all_combinations=True):
    """
    Python version of parser4auto.cpp. 
    Emit all substrings of between minLength and maxLength characters from lines in filenames 
    where the endBasePairs pairs at either end match with each other (Watson-Crick style)
    """

    output_filenames = []
    parsed_sequences = set()
    output_counter = 1

    for filename in filenames:

        try:
            fasta_file = flexibleOpen(filename)
        except (OSError, IOError) as error:
            logging.error("Could not open %s" % filename)
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
            continue

        output_filename = "%s.%s" % (filename, output_counter)
        try:
            with open(output_filename,'w') as output_file:
                for sequence in parsed_sequences:
                    output_file.write("%s\n" % sequence)
        except (IOError, OSError) as error:
            logging.error("Problems trying to create %s: %s" % (output_filename, error))
            return False

        output_filenames.append(output_filename)

    return output_filenames



def runNewcyk(filenames, probabilities_filename):
    output_filenames = []
    for filename in filenames:
        output_filename = Runner.runCommand('newcyk2', filename, [probabilities_filename], 'grm')
        if output_filename:
            output_filenames.append(output_filename)

    # merge the outputs
    output_filename = os.path.join(Runner.OutputDirectory, 'combined.grm')
    command = ['cat'] + output_filenames
    output_file = open(output_filename,'w')
    exit_code = subprocess.call(command, stdout=output_file)
    output_file.close()
    if exit_code != 0:
        logging.error("Problem concatenating together %s into %s" % (', '.join(output_filenames), output_filename))
        return None

    return output_filename

def runCutoffPassscore(filename, score):
    return Runner.runCommand('cutoffpassscore', filename, [str(score)], "cof")



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
    output_filename = Runner.runCommand('RNAfold', filename, ['-C', '--noPS'], 'rfold', output_is_stdout=True, manual_parameters=True, local=False, input_filename=filename)
    return output_filename


def removeMeanFreeEnergyValues(filename):
    """
    Remove the scores produced by RNAfold
    """
    return Runner.runCommand('gawk', filename, ['{print $1}', filename], 'nodg', output_is_stdout=True, manual_parameters=True, local=False)



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


def structuresToDiagrams(filename):
    """
    Draw diagrams for a file of dot-bracket structures 
    """
    output_filename = "%s.diags" % filename

    input_file = open(filename)
    output_file = open(output_filename,'w')

    for line in input_file:
        stripped_line = line.strip()
        if not stripped_line or stripped_line.startswith("#"):
            continue

        if stripped_line.startswith('>'):
            output_file.write(line)
        elif stripped_line[0] in 'ACUGacug':
            sequence = stripped_line
            output_file.write(line)
            structure = next(input_file).strip()
            diagram = diagramFromStructure(sequence, structure)
            output_file.write("%s\n" % "\n".join(diagram))

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
    parser.add_argument("-dG", "--upper-dg", dest="upperDGCutoff", type=float, default=DefaultUpperDGCutoff,
        help="Upper dG cutoff value (default: %(default)s)")
    parser.add_argument("-b", "--end-base-pairs", dest="endBasePairs", type=int, default=DefaultEndBasePairs,
        help="Number of bases to force pair at the ends (default: %(default)s)")
    parser.add_argument("-1", "--one-candidate", dest="oneCandidatePerLine", default=False,
        action="store_true",
        help="Generate only one candidate miRNA per sequence line instead of all possible subsequences")    
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
    parser.add_argument("--probabilities-filename", dest="probabilitiesFilename", default=DefaultProbabilitiesFilename,
        help="Filename where to find the CFG probabilities (default: %(default)s)")
    parser.add_argument('--sge', dest="sge", default=False, action="store_true",
        help="Use Sun Grid Engine")
    parser.add_argument('--queue', dest="sge_queue", default=None,
        help="SGE queue to use")
    parser.add_argument('--sge-logs', dest="sge_logs", default=None,
        help="Redirect SGE log output to this directory")
    parser.add_argument('--max-processes', dest="max_processes", default=None, type=int,
        help="""Maximum number of parallel processes to use. By default, it will be the number of CPU threads if not using SGE, or the number of slots in the queue if using SGE. If using SGE, but not specifying a queue, it is recommended that you pass this parameter""")
    parser.add_argument('-o', '--output-directory', dest="output_directory", default=Runner.OutputDirectory,
        help="""Where to store all the output (default: %(default)s)""")


    parser.add_argument("sequences", nargs="+",
                      help="Filenames containing FASTA or FASTQ sequences")


    parameters = parser.parse_args(args)
    # logger doubles as a section timer
    logging.basicConfig(level=levelFromVerbosity(parameters.verbosity), format="%(asctime)s:%(levelname)s:%(name)s:%(message)s")

    if parameters.sge:
        # Reroute all command-calling to SGE
        import sge
        if parameters.sge_logs:
            sge_logs_directory = parameters.sge_logs
            if not os.path.isdir(sge_logs_directory):
                if os.path.exists(sge_logs_directory):
                    parser.error("%s exists but isn't a directory" % sge_logs_directory)
                else:
                    try:
                        os.makedirs(sge_logs_directory)
                    except OSError as error:
                        parser.error("Problems creating %s: %s" % (sge_logs_directory, error))
        else:
            sge_logs_directory = None

        Runner.Runner = sge.SGE(queue=parameters.sge_queue, logs_directory=sge_logs_directory)

    if parameters.max_processes:
        Runner.MaxProcesses = parameters.max_processes

    Runner.OutputDirectory = os.path.expanduser(parameters.output_directory)


    logging.info("Max processes that will be used: %s" % Runner.max_processes())

    result = 0

    full_filenames = []
    for filename in parameters.sequences:
        full_filename = os.path.expanduser(filename)
        if not os.path.isfile(full_filename):
            parser.error("%s is not a file" % filename)
        full_filenames.append(full_filename)        

    start_time = time.time()

    parsed_filenames = generatePossibleSubsequences(full_filenames, parameters.endBasePairs, parameters.minLength, parameters.maxLength, all_combinations=not parameters.oneCandidatePerLine)

    if not parsed_filenames:
        logging.error("Couldn't parse sequences from %s." % ', '.join(parameters.sequences))
        return 1

    grammar_filename = runNewcyk(parsed_filenames, parameters.probabilitiesFilename)
    if not grammar_filename:
        logging.error("Grammar run failed on %s" % ', '.join(parsed_filenames))
        return 1

    cutoff_filename = runCutoffPassscore(grammar_filename, parameters.grammarCutoff)
    if not cutoff_filename:
        logging.error("Cuttoff run failed on %s" % grammar_filename)
        return 1

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

    diagram_filename = structuresToDiagrams(merged_loops_filename)
    if not diagram_filename:
        logging.error("Problems converting structures to diagrams on %s" % merged_loops_filename)
        return 1


    scores_filename = Runner.runCommand('Scores4mStruct', diagram_filename, [], 'scores')
    if not scores_filename:
        logging.error("Problems running Scores4mStruct on %s" % diagram_filename)
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
