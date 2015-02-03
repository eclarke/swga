# -*- coding: utf-8 -*-
"""primers.py

Functions for parsing primers from a file.

"""
from __future__ import with_statement, division
import subprocess
import csv
import os
import re
import struct
import swga
import swga.resources as resources
import peewee as pw


# The primer database must be initialized before use
# ex: `db.init(db_fname)`
db_fname = 'primer.db'
db = pw.SqliteDatabase(None)


class Primer(pw.Model):
    pid = pw.IntegerField(default=-1)
    seq = pw.TextField(unique=True)
    fg_freq = pw.IntegerField(default=0)
    bg_freq = pw.IntegerField(default=0)
    ratio = pw.FloatField(default=0.0)
    tm = pw.FloatField(default=0.0)
    locations = pw.TextField(default="")
    active = pw.BooleanField(default=False)

    class Meta:
        database = db

    def __repr__(self):
        rep_str = "Primer {0}:{1} (fg_freq:{2}, bg_freq:{3}, ratio:{4})"
        return rep_str.format(
            self.id, self.seq, self.fg_freq, self.bg_freq, self.ratio)


def count_kmers(k, genome_fp, cwd, threshold=1):
    assert isinstance(threshold, int)
    dsk = resources.get_dsk()
    genome = genome_fp.split(os.sep).pop()
    out = '%s-%dmers' % (genome, k)
    outfile = os.path.join(cwd, out + '.solid_kmers_binary')
    if os.path.isfile(outfile):
        swga.message("Binary kmer file already found at %s, skipping..."
                     % outfile)
    else:
        cmdstr = ("{dsk} {genome_fp} {k} -o {out} -t {threshold}"
                  .format(**locals()))
        swga.message("In {cwd}:\n> {cmdstr}".format(**locals()))
        subprocess.check_call(cmdstr, shell=True, cwd=cwd)
    print threshold
    primers = dict((kmer, freq)
                   for kmer, freq in parse_kmer_binary(outfile)
                   if freq >= threshold)
    return primers
        

def parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        kmer_nbits = struct.unpack('i', f.read(4))[0]
        k = struct.unpack('i', f.read(4))[0]
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits // 8),
                                            f.read(kmer_nbits // 8))
                freq = struct.unpack('I', f.read(4))[0]

                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[(kmer_binary[i//4] >> (2 * (i%4))) % 4] + kmer

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


def read_primer_list(plist):
    '''
    Reads in a list of primers, one per line, and returns the corresponding
    records from the primer database.
    '''
    seqs = (re.split(r'[ \t]+', line.strip('\n'))[0] for line in plist)
    return Primer.select().where(Primer.seq << seqs).execute()


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


