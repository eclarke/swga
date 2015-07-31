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
import peewee as pw




def activate(primers):
    '''Marks a list of primers as active.'''
    n = Primer.update(active=True).where(
        Primer.seq << primers).execute()
    return n


def update_locations(primers, fg_genome_fp):
    '''
    Updates the primers from the given set who are missing location data.
    '''
    targets = list(
        Primer.select()
        .where(
            (Primer.seq << primers) &
            (Primer._locations >> None)))
    for primer in targets:
        primer._update_locations(fg_genome_fp)
    swga.database.update_in_chunks(targets, label="Updating primer db... ")


def update_Tms(primers):
    targets = list(
        Primer.select()
        .where(
            (Primer.seq << primers) &
            (Primer.tm >> None)))
    for primer in targets:
            primer.update_tm()
    swga.database.update_in_chunks(targets, label="Updating primer db... ")


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
        cmdstr = ("{dsk} {genome_fp} {k} -o {out} -t {threshold}"
                  .format(**locals()))
        swga.message("In {cwd}:\n> {cmdstr}".format(**locals()))
        try:
            subprocess.check_call(cmdstr, shell=True, cwd=cwd)
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
                    kmer = "ACTG"[(kmer_binary[i//4] >> (2 * (i%4))) % 4] + kmer
                yield kmer, freq
        except struct.error:
            pass


def read_primer_list(lines, fg_genome_fp, bg_genome_fp):
    '''
    Reads in a list of primers, one per line, and returns the corresponding
    records from the primer database. If the primer doesn't exist in the db,
    tries to create it manually. If the primer doesn't appear in the fg genome,
    it skips it with a warning.
    '''
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in lines]
    primers = list(Primer.select().where(Primer.seq << seqs).execute())
    if len(primers) < len(seqs):
        primer_seqs = [p.seq for p in primers]
        missing = [_ for _ in seqs if _ not in primer_seqs]
        for seq in missing:
            swga.message(seq + " not in the database; skipping. Add it "
                         "manually with `swga count --input <file>` ")
    return primers


def parse_kmer_file(lines):
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in lines]
    return seqs


def max_consecutive_binding(mer1, mer2):
    '''
    Return the maximum number of consecutively binding mers
    when comparing two different mers, using the reverse compliment.
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


