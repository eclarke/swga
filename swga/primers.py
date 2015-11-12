# -*- coding: utf-8 -*-
"""primers.py

"""
from __future__ import with_statement, division
import subprocess
import os
import re
import struct
import swga
import swga.database
from swga.database import Primer
import swga.locate
import swga.utils.resources as resources
import swga.melting


def count_kmers(k, genome_fp, cwd, threshold=1):
    assert isinstance(threshold, int)
    dsk = resources.get_dsk()
    genome = genome_fp.split(os.sep).pop()
    out = '%s-%dmers' % (genome, k)
    outfile = os.path.join(cwd, out + '.solid_kmers_binary')
    if os.path.isfile(outfile):
        swga.message("Binary kmer file already exists at %s; parsing..."
                     % outfile)
    else:
        cmdarr = [dsk, genome_fp, str(k), '-o', out, '-t', str(threshold)]
        cmdstr = " ".join(cmdarr)
        swga.message("In {cwd}:\n> {cmdstr}".format(**locals()))
        try:
            subprocess.check_call(cmdarr, cwd=cwd)
        except:
            if os.path.isfile(outfile):
                os.remove(outfile)
            raise
    primers = dict((kmer, freq)
                   for kmer, freq in _parse_kmer_binary(outfile)
                   if freq >= threshold)
    return primers


def _parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        try:
            kmer_nbits = struct.unpack('i', f.read(4))[0]
            k = struct.unpack('i', f.read(4))[0]
        except struct.error:
            if os.path.isfile(fp):
                os.remove(fp)
            raise
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits // 8),
                                            f.read(kmer_nbits // 8))
                freq = struct.unpack('I', f.read(4))[0]
                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[
                        (kmer_binary[i // 4] >> (2 * (i % 4))) % 4] + kmer
                yield kmer, freq
        except struct.error:
            pass


def parse_kmer_file(lines):
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in lines]
    return seqs
