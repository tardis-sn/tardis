#!/usr/bin/env python
from tardis.io import config_reader

from tardis.simulation import Simulation
import numpy as np

import logging
import argparse

import warnings



tardis_description =\
"""
TARDIS Supernova Montecarlo Radiative transfer code

"""

parser = argparse.ArgumentParser(description=tardis_description)
argparse.ArgumentParser()
parser.add_argument('config_fname', help='path to the configuration yaml file')
parser.add_argument('spectrum_fname', help=
'path where to write the output spectrum to. If a virtual spectrum is '
'requested, only the virtual spectrum will be written')


parser.add_argument('--log_file', default=None, help="Name of the log file "
                                                     "(not implemented yet)")

parser.add_argument('--packet_log_file', default=None, help=
"Name of the packet log file. Packet logging needs to be switched on before "
"compiling.")

parser.add_argument('--profile', action='store_true', help=
'run tardis in profiling mode and output results to a specified log file')
parser.add_argument('--profiler_log_file', default='profiler.log', help=
'name of the profiler output file')

parser.add_argument('--gdb', action='store_true', help='print pid and pause')

args = parser.parse_args()

packet_logging_fname = 'tardis_packets.log'

if args.log_file:
    logger = logging.getLogger('tardis')
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

if args.packet_log_file:
    logger = logging.getLogger('tardis_packet_logger')
    logger.setLevel(logging.DEBUG)
    packet_logging_handler = logging.FileHandler(packet_logging_fname, mode='w')
    packet_logging_handler.setLevel(logging.DEBUG)
    packet_logging_formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(packet_logging_formatter)
    logger.addHandler(packet_logging_handler)

tardis_config = config_reader.Configuration.from_yaml(args.config_fname)
simulation = Simulation.from_config(tardis_config)


def get_virtual_spectrum():
    # Catch warning when acessing invalid spectrum_virtual
    with warnings.catch_warnings(record=True) as w:
        spectrum = simulation.runner.spectrum_virtual

        if len(w) > 0 and w[-1]._category_name == 'UserWarning':
            warnings.warn(
                    'Virtual spectrum is not available, using the '
                    'real packet spectrum instead.')
            spectrum = simulation.runner.spectrum
    return spectrum

print('Saving the {} spectrum.'.format(tardis_config.spectrum.method))

# Believe it or not, that is a fancy switch-case statement
# We need functions so that only one spectrum is accessed so we can catch
# the warning properly
get_spectrum = {
        'real': lambda: simulation.runner.spectrum,
        'virtual': get_virtual_spectrum,
        'integrated': lambda: simulation.runner.spectrum_integrated,
        }[tardis_config.spectrum.method]


if args.gdb:
    import os
    print(os.getpid())
    raw_input()  # Workaround to attach gdb
if args.profile:
    import cProfile
    cProfile.runctx('simulation.run()', locals(),
                    globals(), filename=args.profiler_log_file)
else:
    simulation.run()

get_spectrum().to_ascii(args.spectrum_fname)
