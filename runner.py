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
    def __init__(self, *args, **kwargs):
        output = kwargs.pop('output', None)
        error = kwargs.pop('error', None)
        input = kwargs.pop('input', None)
        list.__init__(self, *args, **kwargs)
        self.input = input
        self.error = error
        self.output = output



class LocalRunner(object):
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

        output_filename = output_filename or getattr(command, 'output', None)
        if output_filename:
            output_file = open(output_filename, 'w')
            parameters['stdout'] = output_file

        input_filename = input_filename or getattr(command, 'input', None)
        if input_filename:
            input_file = open(input_filename)
            parameters['stdin'] = input_file

        error_filename = error_filename or getattr(command, 'error', None)
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


    def multi_run(self, commands, environment=None, max_simultaneous=None, buffer_size=16384):
        """
        Run commands and wait for them to finish
        """
        import asyncproc

        DefaultParameters = {
            'close_fds' : True,
            'stderr' : subprocess.PIPE,
            'stdout' : None,
            'stdin' : None
        }

        success = True
        if environment:
            DefaultParameters['env'] = environment

        if buffer_size is not None:
            DefaultParameters['bufsize'] = buffer_size

        processes = {}
        while True:
            finished_in_cycle = False

            while commands and ((max_simultaneous is None) or (len(processes) < max_simultaneous)):
                # we can generate more
                command = commands.pop()
                parameters = DefaultParameters.copy()

                message_parts = ["Running '%s'" % " ".join(command)]

                output_filename = getattr(command, 'output', None)
                input_filename = getattr(command, 'input', None)
                error_filename = getattr(command, 'error', None)
                
                if input_filename:
                    input_file = open(input_filename)
                    parameters['stdin'] = input_file
                    message_parts.append(' input from %s' % input_filename)
                else:
                    input_file = None

                if output_filename:
                    output_file = open(output_filename, 'w')
                    parameters['stdout'] = output_file
                    message_parts.append(' output to %s' % output_filename)
                else:
                    output_file = None

                if error_filename:
                    message_parts.append(' error to %s' % error_filename)
                    if error_filename != output_filename:
                        error_file = open(error_filename, 'w')
                    else:
                        error_file = output_file

                    parameters['stderr'] = error_file
                else:
                    error_file = None

                logging.info(' '.join(message_parts))
                process = asyncproc.Process(command, **parameters)
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


