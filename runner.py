#coding: utf8

from __future__ import print_function, absolute_import, division, unicode_literals

import os, sys
import logging
import subprocess
import time

class Command(list):
    """
    Wafer-thin wrapper around a list so that we can specify an output, an error and an input
    """

    # Marker denoting a boundary between two processes where the preceding process sends the output to the following process
    class PipeMarker(object):
        pass
    PIPE_MARKER = PipeMarker()  

    def __init__(self, *args, **kwargs):
        output = kwargs.pop('output', None)
        error = kwargs.pop('error', None)
        input = kwargs.pop('input', None)
        list.__init__(self, *args, **kwargs)
        self.input = input
        self.error = error
        self.output = output



class Runner(object):
    """
    Base class for objects that take care of running processes
    """

    def split_commands(self, command):
        chain = []
        process_command = []
        for part in command:
            if part is Command.PIPE_MARKER:
                chain.append(process_command)
                process_command = []
            else:
                process_command.append(part)
        if process_command:
            chain.append(process_command)

        return chain


class LocalRunner(Runner):
    """
    Runs processes locally
    """

    DefaultPollingPeriod = 30

    def __init__(self, polling_period=None):
        if polling_period is not None:
            self.polling_period = polling_period
        else:
            self.polling_period = self.DefaultPollingPeriod


    def run(self, command, output_filename=None, error_filename=None, input_filename=None, environment=None):
        parameters = {
            'close_fds' : True
        }

        if environment:
            parameters['env'] = environment

        chain = self.split_commands(command)

        output_filename = output_filename or getattr(command, 'output', None)
        if output_filename:
            output_file = open(output_filename, 'w')
        else:
            output_file = None

        input_filename = input_filename or getattr(command, 'input', None)
        if input_filename:
            input_file = open(input_filename)
        else:
            input_file = None

        error_filename = error_filename or getattr(command, 'error', None)
        if error_filename:
            if error_filename != output_filename:
                error_file = open(error_filename, 'w')
            else:
                error_file = output_file
        else:
            error_file = None


        last_out = input_file
        intermediate_processes = []
        for intermediate_command in chain[:-1]:
            parameters['stdin'] = last_out
            parameters['stderr'] = None
            parameters['stdout'] = subprocess.PIPE

            logging.debug("Running %s" % " ".join(intermediate_command))
            process = subprocess.Popen(intermediate_command, **parameters)
            intermediate_processes.append(process)
            last_out = process.stdout

        output_command = chain[-1]
        parameters['stdin'] = last_out
        parameters['stdout'] = output_file
        parameters['stderr'] = error_file

        logging.debug("Running %s" % " ".join(output_command))
        exit_code = subprocess.call(output_command, **parameters)
        if exit_code == 0:
            logger = logging.debug
        else:
            logger = logging.error
        logger("Finished running %s with exit code %s" % (" ".join(output_command), exit_code))

        if input_file is not None:
            input_file.close()            
        if output_file is not None:
            output_file.close()
        if error_file is not None and not error_file.closed:
            error_file.close()

        for process in intermediate_processes:
            process.wait()

        return exit_code


    def multi_run(self, commands, environment=None, max_simultaneous=None, buffer_size=16384):
        """
        Run commands and wait for them to finish
        """
        import asyncproc

        OutputParameters = {
            'close_fds' : True,
            'stderr' : subprocess.PIPE,
            'stdout' : None,
            'stdin' : None
        }

        InputParameters = {
            'close_fds' : True,
            'stderr' : None,
            'stdout' : subprocess.PIPE,
            'stdin' : None
        }


        success = True
        if environment:
            InputParameters['env'] = environment
            OutputParameters['env'] = environment

        if buffer_size is not None:
            InputParameters['bufsize'] = buffer_size
            OutputParameters['bufsize'] = buffer_size



        processes = {}
        while True:
            finished_in_cycle = False

            while commands and ((max_simultaneous is None) or (len(processes) < max_simultaneous)):
                # we can generate more
                command = commands.pop()
                chain = self.split_commands(command)

                file_message_parts = []

                output_filename = getattr(command, 'output', None)
                input_filename = getattr(command, 'input', None)
                error_filename = getattr(command, 'error', None)
                
                if input_filename:
                    input_file = open(input_filename)
                    file_message_parts.append(' input from %s' % input_filename)
                else:
                    input_file = None

                if output_filename:
                    output_file = open(output_filename, 'w')
                    file_message_parts.append(' output to %s' % output_filename)
                else:
                    output_file = None

                if error_filename:
                    file_message_parts.append(' error to %s' % error_filename)
                    if error_filename != output_filename:
                        error_file = open(error_filename, 'w')
                    else:
                        error_file = output_file
                else:
                    error_file = None

                last_input = input_file

                message_parts = []
                for process_command in chain[:-1]:
                    # everything but the output

                    parameters = InputParameters.copy()
                    message_parts.append("'%s'" % process_command)

                    parameters['stdin'] = last_input

                    process = subprocess.Popen(process_command, **parameters)
                    last_input = process.stdout

                output_command = chain[-1]
                message_parts.append("'%s'" % output_command)

                parameters = OutputParameters.copy()
                parameters['stdin'] = last_input
                parameters['stdout'] = output_file
                if error_file is not None:
                    parameters['stderr'] = error_file

                logging.info("Running %s %s" % (' '.join(message_parts), ' '.join(file_message_parts)))
                process = asyncproc.Process(output_command, **parameters)
                processes[process] = (command, input_file, output_file, error_file)

            for process, command_details in processes.items():
                if command_details[-1] is None:
                    # if we are not piping stderr to a file, read it and log it
                    error = process.readerr()
                    if error:
                        logging.info("Process %s err: %s" % (process.pid(), error))

                exit_code = process.wait(os.WNOHANG) # check if we are done
                if exit_code is not None:
                    command, input_file, output_file, error_file = command_details
                    if input_file is not None:
                        input_file.close()

                    if output_file is not None:
                        output_file.close()

                    if error_file is not None and not error_file.closed:
                        error_file.close()

                    finished_in_cycle = True

                    del processes[process]
                    if exit_code != 0:
                        logger = logging.error
                        success = False
                    else:
                        logger = logging.info
                    logger("'%s' exited with code %s" % (' '.join(command), exit_code))

            # if there's nothing left, quit

            if not commands and not processes:
                break

            # if nothing finished, sleep for a while, otherwise go for more
            if not finished_in_cycle:
                time.sleep(self.polling_period)


        return success


    def max_processes(self):
        import multiprocessing
        return multiprocessing.cpu_count()


