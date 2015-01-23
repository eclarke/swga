# -*- coding: utf-8 -*-
"""primers.py

Functions for parsing primers from a file.

"""
import re
import sys
from collections import namedtuple
from .core import errprint

Primer = namedtuple('Primer', 'id, seq, bg_freq, fg_freq, ratio')

def parse_primer(string, line_no=1):
    '''
    Takes a line from a tab- or space-delimited file where each row specifies
    a primer sequence, fg binding count, bg binding count, and fg/bg binding
    count ratio.

    Returns: A Primer object (or raises ValueError if it cannot parse the line)
    '''
    seq, fg_freq, bg_freq, ratio = re.split(r'[ \t]+', string.strip('\n'))
    return Primer(line_no, seq, int(bg_freq), int(fg_freq), float(ratio))


def read_primer_file(file_handle, echo_input=False):
    '''
    Calls parse_primer() on each line of the input file. If a malformed line is
    found, will skip parsing and (optionally) output a warning message.

    Arguments:
    file_handle: an open file handle to read from
    echo_input:  if True, echo each line of the file to stdout
    verbose:     if True, warn when parsing a line fails

    Returns: A list of Primer objects
    '''
    primers = []
    for i, line in enumerate(file_handle):
        try:
            primers.append(parse_primer(line, i+1))
            if echo_input:
                sys.stdout.write(line)
        except ValueError as e:
            if verbose:
                errprint("Skipping line %i (reason: %s)\n" % (i, e.message))
            continue
    return primers
