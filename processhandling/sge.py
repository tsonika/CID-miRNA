"""
Module that handles the interaction with the Sun Grid Engine to run programs on a cluster
with shared drive space
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, sys, logging
import subprocess
import time
import re

from processhandling.runner import Runner, Command

_find_unsafe = re.compile(r'[^\w@%+=:,./-]').search

def quote(s):
    """
    (Taken from Python 3's shlex module)
    Return a shell-escaped version of the string *s*.
    """
    if not s:
        return "''"
    if _find_unsafe(s) is None:
        return s

    # use single quotes, and put single quotes into double quotes
    # the string $'b is then quoted as '$'"'"'b'
    return "'" + s.replace("'", "'\"'\"'") + "'"



class SGE(Runner):
    MaxProcessesNoQueue = 10
    DefaultPollingPeriod = 30

    def __init__(self, queue=None, logs_directory=None, polling_period=DefaultPollingPeriod):
        """
        If logs_directory is given, stdout and stderr will be redirected there if nothing else
        is given
        """

        self.queue = queue
        self.polling_period = polling_period
        self.logs_directory = logs_directory

    def run(self, command, queue=None, name=None, sync=True, output_filename=None, error_filename=None, input_filename=None, environment=None, append=False):
        """
        Run a command on given queue. If sync is True, don't return until the command is done
        Returns command's exit code if in sync, SGE's process ID otherwise, or -1 if the SGE call
        failed
        """

        output_filename = output_filename or getattr(command, 'output', None)
        input_filename = input_filename or getattr(command, 'input', None)
        error_filename = error_filename or getattr(command, 'error', None)

        sge_command = self.wrap(command, sync=sync, queue=queue, name=name, output_filename=output_filename,
            error_filename=error_filename, input_filename=input_filename)

        if not append:
            if output_filename and os.path.exists(output_filename):
                try:
                    os.remove(output_filename)
                except OSError as error:
                    logging.error("Problems removing %s: %s" % (output_filename, error))
                    return -1
            if error_filename and os.path.exists(error_filename):
                try:
                    os.remove(error_filename)
                except OSError as error:
                    logging.error("Problems removing %s: %s" % (error_filename, error))

        parameters = {
            'close_fds' : True
        }
        if environment:
            parameters['env'] = environment

        if sync:
            logging.info("Calling %s" % ' '.join(sge_command))
            exit_code = subprocess.call(sge_command, **parameters)
            logging.info("Finished running %s with exit code %s" % (" ".join(sge_command), exit_code))
            return exit_code
        else:
            try:
                logging.info("Calling (async) %s" % ' '.join(sge_command))
                output = subprocess.check_output(sge_command, **parameters)
            except subprocess.CalledProcessError:
                return -1
            else:
                try:
                    job_id = int(output)
                except ValueError:
                    logging.error("Submitted %s to SGE, but didn't get a job ID. Got %s" % (' '.join(command), output))
                    return -1
                else:
                    logging.info("SGE process ID: %s" % job_id)
                    return job_id


    def multi_run(self, commands, queue=None, polling_period=None, environment=None, max_simultaneous=None):
        """
        Run and wait for commands to finish. Return if they all succeeded
        commands should be a list of commands, each of which is an iterable
        """

        if polling_period is None:
            polling_period = self.polling_period

        if queue is None:
            queue = self.queue

        processes = {}

        allOK = True
        for command in commands:
            sge_process_id = self.run(command, queue=queue, sync=False, environment=environment)
            if sge_process_id < 0:
                allOK = False
                break
            processes[sge_process_id] = command
        
        if allOK:
            while True:
                for process_id, command in processes.items():
                    exit_code = self.get_exit_code(process_id)
                    if exit_code is not None:
                        if exit_code != 0:
                            logger = logging.error
                            allOK = False
                        else:
                            logger = logging.info

                        logger("'%s' exited with code %s" % (' '.join(argument if argument is not Command.PIPE_MARKER else '|' for argument in command), exit_code))

                        del processes[process_id]
                if not processes:
                    break
                time.sleep(polling_period)

        return allOK


    def get_exit_code(self, job_id):
        """
        Get the exit code for an SGE job id
        """

        command = ['qacct', '-j', str(job_id)]
        try:
            output = subprocess.check_output(command, stderr=subprocess.STDOUT, close_fds=True)
        except subprocess.CalledProcessError:
            # process didn't exist
            return None
        else:
            for line in output.splitlines():
                if line.startswith('exit_status'):
                    parts = line.split()
                    try:
                        exit_code = int(parts[-1].strip())
                    except ValueError:
                        logging.error("Looking for exit code for SGE process %s, got %s" % (job_id, line))
                        return None
                    else:
                        return exit_code

    def wrap(self, command, sync=True, name=None, queue=None, output_filename=None, error_filename=None,
        input_filename=None):
        """
        Wrap a command so that it can run on SGE
        """

        if queue is None:
            queue = self.queue

        sge_command = 'qsub -terse -cwd -b Y -V -S /bin/sh'.split()
        if sync:
            sge_command.extend(['-sync', 'yes'])
        if name:
            sge_command.extend(['-N', name])
        if queue:
            sge_command.extend(['-q', queue])
        if output_filename:
            sge_command.extend(['-o', output_filename])
        elif self.logs_directory:
            sge_command.extend(['-o', self.logs_directory])
        if error_filename:
            if error_filename != output_filename:
                sge_command.extend(['-e',error_filename])
            else:
                sge_command.extend(['-j', 'yes'])
        elif self.logs_directory:
            sge_command.extend(['-e', self.logs_directory])
        if input_filename:
            sge_command.extend(['-i', input_filename])


        chain = self.split_commands(command)
        # manual list joined by '|'
        for process_command in chain[:-1]:
            sge_command.extend([quote(argument) for argument in process_command])
            sge_command.append('|')
        sge_command.extend([quote(argument) for argument in chain[-1]])
        return sge_command


    def max_processes(self, queue=None):
        """
        Return the number of slots in the queue
        """
        if queue is None:
            queue = self.queue

        if not queue:
            # It's hard to find which queue will be used by default, so just go with some
            # default number
            return self.MaxProcessesNoQueue

        command = ['qconf', '-sq', queue]
        try:
            output = subprocess.check_output(command, close_fds=True)
        except subprocess.CalledProcessError:
            logging.error("Couldn't get maximum number of processes")
            return None

        for line in output.splitlines():
            if line.startswith('slots'):
                bits = line.split()
                try:
                    slots_number = bits[1].split(',')[0]
                    slots = int(slots_number)
                except (ValueError, IndexError):
                    logging.error("Slots line didn't come in the format we expected: %s" % line)
                    return None

                return slots
        else:
            logging.error("Couldn't find the number of slots for the SGE queue")
            return None
