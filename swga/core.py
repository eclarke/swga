# -*- coding: utf-8 -*-
"""core.py

This module contains utility functions and constructors for other modules in
SWGA code. In particular, it contains the default config parser and command-line
argument parsers for the subcommands, as well as error and warning functions.

"""

import os
import sys
import json
import errno
import argparse
import textwrap
import ConfigParser

from clint.textui import puts, colored, STDERR, indent

default_config_file = 'parameters.cfg'



# def get_swgahome():
#     '''
#     Gets the SWGAHOME environmental variable, catching common
#     errors while doing so.
#     '''
#     swgahome = os.environ.get('SWGAHOME')
#     if not swgahome:
#         swga_error("SWGAHOME environment variable not set.")
#     if not os.path.isabs(swgahome):
#         swga_error("SWGAHOME must be an absolute path.")
#     return swgahome


def mkdirp(path):
    '''Simulates 'mkdir -p': creates a directory unless it already exists'''
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def parse_config(cfg_file, section):
    '''
    Parses a config file and returns a dictionary of the values found
    in the specified section, along with the ConfigParser itself.
    '''
    config = ConfigParser.SafeConfigParser()
    defaults = {}
    with open(cfg_file) as cfg_file_fp:
        config.readfp(cfg_file_fp)
        defaults = dict(config.items(section))
        if not all(defaults.values()):
            for key, value in defaults.iteritems():
                if not value:
                    swga_warn("""
                    Warning: value for `{0}` unspecified in config file.""".format(key))
        return defaults, config


def basic_cmd_parser(description, cmd_name, cfg_file):
    try:
        defaults, _ = parse_config(cfg_file, cmd_name)
    except IOError:
        defaults = dict()
    parser = argparse.ArgumentParser(description=description, prog='swga '+cmd_name)
    parser.set_defaults(**defaults)
    return parser

def swga_error(msg, errcode=1):
    '''Prints an error message to stderr and exits.'''
    errprint('{}\n'.format(msg))
    sys.exit(errcode)


def swga_warn(msg):
    '''Prints a warning message to stderr.'''
    with indent(2, quote=colored.red("!! ")):
        puts(textwrap.fill(textwrap.dedent(msg)), stream=STDERR)


def errprint(text):
    puts(colored.red(textwrap.fill(textwrap.dedent(text)).strip()),
         stream=STDERR)


def print_status(prog_name, args, cfg_file, from_stdin):
    
    puts("Command: {}".format(prog_name), stream=STDERR)
    with indent(4, quote="# "):
        puts("Config file: {}".format(os.path.abspath(cfg_file)),
             stream=STDERR)
        puts("Parameters: \n{}".format(json.dumps(vars(args), 
                                                  sort_keys=True, 
                                                  indent=2,
                                                  separators=(',', ': '),
                                                  default=lambda x: x.name)),
                 stream=STDERR)
        if from_stdin:
            puts("Receiving input from stdin...")


def progressbar(i, length):
    if i >= 1:
        i = i/(length*1.0)
    sys.stderr.write('\r[%-20s] %-3d%%' % ('='*int(round(i*20)), i*100))
    sys.stderr.flush()
