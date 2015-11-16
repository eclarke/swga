# -*- coding: utf-8 -*-
'''Functions for counting/reading raw k-mer sequences.'''
from __future__ import with_statement, division
import subprocess
import os
import re
import struct

from swga import message
#from swga.utils import dsk
import swga.utils


def parse_kmer_file(lines):
    '''Returns the first field of whitespace-delimited lines.

    :param lines: an iterable of strings (e.g. an open handle)
    '''
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in lines]
    return seqs


def count_kmers(k, genome_fp, cwd, threshold=1):
    '''Counts k-mers in the specified genome.

    :param k: the number of nucleotides in the k-mers
    :param genome_fp: the file path to the genome/fasta file
    :param cwd: the current working directory (to store intermediate cache)
    :param threshold: the minimum k-mer frequency
    '''

    assert isinstance(threshold, int)

    genome = genome_fp.split(os.sep).pop()
    out = '%s-%dmers' % (genome, k)
    outfile = os.path.join(cwd, out + '.solid_kmers_binary')

    if os.path.isfile(outfile):
        message("Binary kmer file already exists at %s; parsing..."
                % outfile)
    else:
        cmdarr = [swga.utils.dsk, genome_fp, str(k), '-o', out, '-t', str(threshold)]
        cmdstr = " ".join(cmdarr)
        message("In {cwd}:\n> {cmdstr}".format(**locals()))
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


def max_sequential_nt(mer1, mer2):
    '''Returns the maximum number of nt that would bind between two kmers.

    :param kmer1: the first kmer
    :param kmer2: the second kmer (order unimportant)
    '''

    binding = { 'A': 'T', 'T': 'A',
                'C': 'G', 'G': 'C',
                '_':  False}

    # Swap variables if the second is longer than the first
    if len(mer2) > len(mer1):
        mer1, mer2 = mer2, mer1

    # save the len because it'll change when we do a ljust
    mer1_len = len(mer1)
    # reverse mer2,
    mer2 = mer2[::-1]
    # pad mer one to avoid errors
    mer1 = mer1.ljust(mer1_len + len(mer2), "_")

    max_bind = 0
    for offset in range(mer1_len):
        consecutive = 0
        for x in range(len(mer2)):
            if binding[mer1[offset+x]] == mer2[x]:
                consecutive += 1
                if consecutive > max_bind:
                    max_bind = consecutive
            else:
                consecutive = 0
    return max_bind


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
