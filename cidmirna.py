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
from foldingsubsequences import DefaultMinLength, DefaultMaxLength, DefaultEndBasePairs
from utils import generateRNACombinations
from runner import LocalRunner, Command

DefaultUpperDGCutoff = -13.559
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


class Configuration(object):
    """
    Class/Singleton that makes it easier to configure what method is used to run external programs
    and where the output should go
    """

    Runner = LocalRunner()
    MaxProcesses = None
    OutputDirectory = '.'

    @classmethod
    def get_command_and_environment(cls, command, local=True):
        ld_library_path = None
        if local:
            script_path = os.path.dirname(__file__)
            if not script_path:
                script_path = '.'

            script_path = os.path.join(script_path, 'bin')

            if command in NeedPreload:
                ld_library_path = script_path
            full_command = [os.path.join(script_path, command)]
            if command.endswith('.py'):
                # prefix the python executable that we are using
                full_command.insert(0, sys.executable)
        else:
            full_command = [command]

        environment = dict(os.environ)
        if ld_library_path:
            environment['LD_LIBRARY_PATH'] = "%s%s" % (ld_library_path, ':%s' % environment.get('LD_LIBRARY_PATH') if "LD_LIBRARY_PATH" in environment else '')

        return full_command, environment

    @classmethod
    def runCommand(cls, command, filename, parameters, output_extension, output_is_stdout=False, input_filename=None, manual_parameters=False, local=True):

        """
        Standard command run wrapper.
        """

        command, environment = cls.get_command_and_environment(command, local=local)

        output_filename = "%s.%s" % (cls.map_file_to_output_directory(filename), output_extension)
        if manual_parameters:
            full_command = command + parameters
        else:
            full_command = command + [filename] + parameters

        if output_is_stdout:
            piped_output = output_filename
        else:
            full_command.append(output_filename)
            piped_output = None

        exit_code = cls.Runner.run(full_command, output_filename=piped_output, 
            input_filename=input_filename, environment=environment)
        if exit_code != 0:
            return False

        return output_filename        


    @classmethod
    def multi_run(cls, commands, local=True):
        for command in commands:
            # patch command if it needs patching
            binary, environment = cls.get_command_and_environment(command[0], local=local)
            command[0:1] = binary

        return cls.Runner.multi_run(commands, environment=environment, max_simultaneous=cls.MaxProcesses)



    @classmethod
    def map_file_to_output_directory(cls, filename):
        return os.path.join(cls.OutputDirectory, os.path.basename(filename))


    @classmethod
    def max_processes(cls):
        if cls.MaxProcesses:
            return cls.MaxProcesses
        else:
            return cls.Runner.max_processes()



def common_prefix(filenames):
    """
    Find a common prefix for a whole bunch of files. 
    """

    basenames = [os.path.basename(filename) for filename in filenames]
    prefix = []
    candidate = basenames[0]
    for index, letter in enumerate(candidate):
        if all(letter == other[index] for other in basenames):
            prefix.append(letter)
        else:
            break

    return ''.join(prefix)


def generatePossibleSubsequencesWrapper(filenames, end_base_pairs, min_length, max_length, all_combinations=True):
    """
    Python version of parser4auto.cpp. 
    Emit all substrings of between minLength and maxLength characters from lines in filenames 
    where the endBasePairs pairs at either end match with each other (Watson-Crick style)
    """

    output_filenames = []
    output_counter = 1

    number_processes = Configuration.max_processes()
    number_files = len(filenames)

    # we want at least two files per process to smooth out the load, but we don't want to 
    # generate more than, say 2000 files

    total_files = float(number_files)
    split_level = 0

    while total_files / number_processes < 2 and total_files < 2000:
        split_level += 1
        total_files *= 4


    logging.info("We are going to split the files to the %dth-level" % split_level)

    if not all_combinations:
        one_only = ['-1']
    else:
        one_only = []

    rnaSequences = generateRNACombinations(split_level)
    commands = []
    for filename in filenames:
        output_filename = "%s.%s" % (Configuration.map_file_to_output_directory(filename), output_counter)

        command = ['findsubsequences.py', '-b', str(end_base_pairs), '-m', str(min_length),
        '-x', str(max_length), '-o', output_filename, '-s', str(split_level)] + one_only + [filename]

        commands.append(command)
        if split_level == 0:
            output_filenames.append(output_filename)
        else:
            output_filenames.extend("%s.%s" % (output_filename, suffix) for suffix in rnaSequences)

    if not Configuration.multi_run(commands):
        logging.error("Some part of the subsequence generation failed")
        return False
    

    return output_filenames


def runUnique(filenames):
    output_filenames = []
    commands = []
    for filename in filenames:
        output_filename = '%s.uniq' % Configuration.map_file_to_output_directory(filename)
        command = Command(['sort', filename, Command.PIPE_MARKER, 'uniq'], output=output_filename)
        output_filenames.append(output_filename)
        commands.append(command)

    if not Configuration.multi_run(commands, local=False):
        logging.error("Some problem with uniquifying")
        return False

    return output_filenames



def runNewcyk(filenames, probabilities_filename):
    output_filenames = []
    commands = []
    for filename in filenames:
        output_filename = '%s.grm' % Configuration.map_file_to_output_directory(filename)
        command = ['newcyk2', filename, probabilities_filename, output_filename]
        output_filenames.append(output_filename)
        commands.append(command)

    if not Configuration.multi_run(commands, local=True):
        logging.error("Some problem with newcyk2")
        return False

    if len(output_filenames) > 1:
        # merge the outputs

        # try to get a reasonable name for the combined file

        prefix = common_prefix(output_filenames)

        if prefix and not prefix.endswith('.'):
            prefix += '.'

        output_filename = Configuration.map_file_to_output_directory('%scombined.grm' % prefix)
        command = ['cat'] + output_filenames
        logging.info("Catting %s into %s" % (', '.join(output_filenames), output_filename))
        output_file = open(output_filename,'w')
        exit_code = subprocess.call(command, stdout=output_file)
        output_file.close()
        logging.info("Finished catting %s into %s" % (', '.join(output_filenames), output_filename))
        if exit_code != 0:
            logging.error("Problem concatenating together %s into %s" % (', '.join(output_filenames), output_filename))
            return None
    else:
        output_filename = output_filenames[0]

    return output_filename

def runCutoffPassscore(filename, score):
    return Configuration.runCommand('cutoffpassscore', filename, [str(score)], "cof")



def convertToFasta(filename):
    output_filename = "%s.fasta" % Configuration.map_file_to_output_directory(filename)

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

    output_filename = "%s.4rfold" % Configuration.map_file_to_output_directory(filename)

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
    output_filename = Configuration.runCommand('RNAfold', filename, ['-C', '--noPS'], 'rfold', output_is_stdout=True, manual_parameters=True, local=False, input_filename=filename)
    return output_filename


def removeMeanFreeEnergyValues(filename):
    """
    Remove the scores produced by RNAfold
    """
    return Configuration.runCommand('gawk', filename, ['{print $1}', filename], 'nodg', output_is_stdout=True, manual_parameters=True, local=False)



def keepOneLoop(structure):
    """
    Keep only the centre loop in a structure. 
    """

    loop_centre = re.compile(r'\(\.*\)')

    
    loop_centres = []
    halfway = len(structure) // 2
    for match in loop_centre.finditer(structure):
        loop_centre = (match.start() + match.end()) // 2
        loop_centres.append((abs(halfway - loop_centre), loop_centre))

    closest_to_centre = min(loop_centres)
    centre_loop = closest_to_centre[1]

    # We've decided on the winner, now get rid of every other loop

    def straightenLoopsBefore(structure, end, closes_loop, opens_loop):
        fixed_until = end
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
                        fixed_until = index
                        break
                index -= 1
        return structure

    # First backwards
    structure = straightenLoopsBefore(structure, centre_loop, ')', '(')
    # Now forwards, by reversing the structure
    structure = structure[::-1]
    structure = straightenLoopsBefore(structure, len(structure) - centre_loop - 1, '(', ')')
    structure = structure[::-1]

    return structure


def mergeLoops(filename):
    """
    Make sure there's only one loop per structure. If there are more, keep the centre one and 
    'straighten' the other ones out
    """

    output_filename = "%s.mloops" % Configuration.map_file_to_output_directory(filename)

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
    output_filename = "%s.diags" % Configuration.map_file_to_output_directory(filename)

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

    output_filename = "%s.pass" % Configuration.map_file_to_output_directory(filename)

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
    output_filename = "%s.fasta" % Configuration.map_file_to_output_directory(filename)

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
    parser.add_argument('-o', '--output-directory', dest="output_directory", default=Configuration.OutputDirectory,
        help="""Where to store all the output (default: %(default)s)""")
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

        Configuration.Runner = sge.SGE(queue=parameters.sge_queue, logs_directory=sge_logs_directory)

    if parameters.max_processes:
        Configuration.MaxProcesses = parameters.max_processes

    Configuration.OutputDirectory = os.path.expanduser(parameters.output_directory)

    if not os.path.exists(Configuration.OutputDirectory):
        logging.info("Creating %s" % Configuration.OutputDirectory)
        os.makedirs(Configuration.OutputDirectory)

    logging.info("Max processes that will be used: %s" % Configuration.max_processes())

    result = 0

    full_filenames = []
    for filename in parameters.sequences:
        full_filename = os.path.expanduser(filename)
        if not os.path.isfile(full_filename):
            parser.error("%s is not a file" % filename)
        full_filenames.append(full_filename)        

    start_time = time.time()

    parsed_filenames = generatePossibleSubsequencesWrapper(full_filenames, parameters.endBasePairs, parameters.minLength, parameters.maxLength, all_combinations=not parameters.oneCandidatePerLine)

    if not parsed_filenames:
        logging.error("Couldn't parse sequences from %s." % ', '.join(parameters.sequences))
        return 1

    uniqued_filenames = runUnique(parsed_filenames)
    if not uniqued_filenames:
        logging.error("Couldn't make unique versions of %s" % ', '.join(parsed_filenames))
        return 1

    grammar_filename = runNewcyk(uniqued_filenames, parameters.probabilitiesFilename)
    if not grammar_filename:
        logging.error("Grammar run failed on %s" % ', '.join(uniqued_filenames))
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


    scores_filename = Configuration.runCommand('Scores4mStruct', diagram_filename, [], 'scores')
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


    if len(full_filenames) > 1:
        common_name = common_prefix(full_filenames)
        if common_name and not common_name.endswith('.'):
            common_name += '.'

        if not common_name:
            common_name = 'combined.'
    else:
        common_name = full_filenames[0]

    final_fasta_filename = Configuration.map_file_to_output_directory("%sfinal.fasta" % common_name)
    os.rename(fasta_filename, final_fasta_filename)

    end_time = time.time()

    logging.info("Took %s seconds" % (end_time - start_time))

    print("FINAL OUTPUT FILES:\n\tPost-grammar cof: %s\n\tPost-structure cof: %s\n" % (grammar_fasta_filename, final_fasta_filename))


    return result

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
