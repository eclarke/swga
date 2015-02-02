# -*- coding: utf-8 -*-
"""primers.py

Functions for parsing primers from a file.

"""
from __future__ import with_statement, division
import subprocess
import csv
import os
import struct
import swga
import swga.resources as resources

class Primer:

    def __init__(self, id, seq, bg_freq=None, fg_freq=None, extra=None):
        self.id = id
        self.seq = seq
        self.bg_freq = int(bg_freq) if bg_freq else 0
        self.fg_freq = int(fg_freq) if fg_freq else 0
        self.ratio = self.fg_freq / self.bg_freq if self.bg_freq is not 0 else 0
        self.extra = extra

    def to_line(self, newline=True):
        primer_str = "{seq}\t{fg_freq}\t{bg_freq}\t{ratio}"
        line = primer_str.format(id=self.id,
                                 seq=self.seq,
                                 fg_freq=self.fg_freq,
                                 bg_freq=self.bg_freq,
                                 ratio=self.ratio)
        line = line + '\n' if newline else line
        return line

    def to_tuple(self):
        return (self.seq, self.fg_freq, self.bg_freq, self.ratio)

    def __repr__(self):
        rep_str = "Primer {0}:{1} (fg_freq:{2}, bg_freq:{3}, ratio:{4})"
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


def count_kmers(k, genome_fp, cwd, threshold=1):
    dsk = resources.get_dsk()
    genome = genome_fp.split(os.sep).pop()
    out = os.path.join(cwd, genome + '.kmers')
    outfile = out + '.solid_kmers_binary'
    cmdstr = "{dsk} {genome} {k} -o {out} -t {threshold}".format(**locals())

    # Run DSK
    swga.mkdirp(cwd)
    subprocess.check_call(cmdstr, shell=True, cwd=cwd)

    return outfile
        

def parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        kmer_nbits = struct.unpack('i', f.read(4))[0]
        k = struct.unpack('i', f.read(4))[0]
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits / 8),
                                            f.read(kmer_nbits / 8))
                freq = struct.unpack('I', f.read(4))[0]
                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[(kmer_binary[i/4] >> (2 * (i%4))) % 4] + kmer
                yield kmer, freq
        except struct.error:
            pass



def read_primer_file(infile, echo_input=False):
    '''
    Calls parse_primer() on each line of the input file. If a malformed line is
    found, will skip parsing and (optionally) output a warning message.

    Arguments:
    infile: an open file handle to read from

    Returns: A list of Primer objects
    '''
    rows = csv.DictReader((f for f in infile if not f.startswith('#')),
                          delimiter='\t',
                          fieldnames=['seq', 'fg_freq', 'bg_freq', 'ratio'], 
                          restkey='extra')
    primers = [Primer(id=i+1, **row) for i, row in enumerate(rows)]
    return primers


def write_primer_file(primers, out_handle, header=True):
    '''Writes each primer in primers to a line.

    Arguments:
    - primers: a list of primers
    - out_handle: an open file handle to write to
    - header: write a header?
    '''
    fields = ["seq", "fg_freq", "bg_freq", "ratio"]
    writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter='\t', 
                            lineterminator='\n',
                            extrasaction='ignore')
    if header:
        out_handle.write("#")
        writer.writeheader()
    [writer.writerow(p.__dict__) for p in primers]


def mk_primer_tbl(connection):
    with connection as c:
        c.execute('''
            create table primers 
            (seq text primary key, fg_freq integer, bg_freq integer, ratio real)
            ''')


def update_primer_tbl(primers, connection):
    primer_tuples = (p.to_tuple() for p in primers)
    with connection as c:
        c.executemany("""
        insert or replace into primers(seq, fg_freq, bg_freq, ratio) 
        values (?,?,?,?)""", primer_tuples)


