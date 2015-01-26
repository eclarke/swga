# -*- coding: utf-8 -*-
"""primers.py

Functions for parsing primers from a file.

"""
import re
import sys
from .core import errprint

class Primer:

    def __init__(self, id, seq, bg_freq=None, fg_freq=None, ratio=None, 
                 old_id=None):
        self.id = id
        self.seq = seq
        self.bg_freq = bg_freq 
        self.fg_freq = fg_freq 
        self.ratio = ratio
        self.old_id = old_id


    def to_line(self, with_id=True, newline=True):
        primer_str = "{seq}\t{fg_freq}\t{bg_freq}\t{ratio}\t{id}"
        line = primer_str.format(id = self.id,
                                 seq=self.seq,
                                 fg_freq = self.fg_freq,
                                 bg_freq = self.bg_freq,
                                 ratio = self.ratio)
        line = line + '\n' if newline else line
        return line


    def __repr__(self):
        rep_str = "Primer #{0}:{1} (fg_freq:{2}, bg_freq:{3}, ratio:{4})"
        return rep_str.format(
            self.id, self.seq, self.fg_freq, self.bg_freq, self.ratio)


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.seq == other.seq and 
                    self.fg_freq == other.fg_freq and
                    self.bg_freq == other.bg_freq and
                    self.ratio == other.ratio)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


def _parse_primer(string, line_no=1):
    '''
    Takes a line from a tab- or space-delimited file where each row specifies
    a primer sequence, fg binding count, bg binding count, and fg/bg binding
    count ratio.

    Returns: A Primer object (or raises ValueError if it cannot parse the line)
    '''
    try:
        # If ID included in line
        seq, fg_freq, bg_freq, ratio, old_id = re.split(r'[ \t]+', string.strip('\n'))
        return Primer(line_no, seq, int(bg_freq), int(fg_freq), float(ratio), old_id)
    except ValueError:
        seq, fg_freq, bg_freq, ratio = re.split(r'[ \t]+', string.strip('\n'))
        return Primer(line_no, seq, int(bg_freq), int(fg_freq), float(ratio))


def read_primer_file(file_handle, echo_input=False):
    '''
    Calls parse_primer() on each line of the input file. If a malformed line is
    found, will skip parsing and (optionally) output a warning message.

    Arguments:
    file_handle: an open file handle to read from
    echo_input:  if True, echo each line of the file to stdout

    Returns: A list of Primer objects
    '''
    primers = []
    for i, line in enumerate(file_handle):
        print i, line
        try:
            primers.append(_parse_primer(line, i+1))
            if echo_input:
                sys.stdout.write(line)
        except ValueError as e:
            errprint("Skipping line %i (reason: %s)\n" % (i, e.message))
            continue
    return primers


def write_primer_file(primers, out_handle, header=True):
    '''Writes each primer in primers to a line.

    Arguments:
    - primers: a list of primers
    - out_handle: an open file handle to write to
    - header: write a header?
    '''
    if header:
        header_str = "#"+"\t".join(["seq", "fg_freq", "bg_freq", "ratio", "id"])
        out_handle.write(header_str+'\n')
    for primer in primers:
        out_handle.write(primer.to_line())
