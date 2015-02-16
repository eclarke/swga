# -*- coding: utf-8 -*-
"""locate.py

This module contains functions for locating primers in genomes (i.e. finding
small substrings in larger strings).

"""
import json
import multiprocessing
import signal
import time

from pyfaidx import Fasta

import swga
from .core import progressbar


def linearize_binding_sites(primers, chr_ends):
    '''
    Modifies the primer binding site locations as if they were positions on a
    linear genome composed of all the chromosomes concatenated together. This
    allows us to compute distances between primer binding sites correctly.
    '''
    new_locs = []
    for primer in primers:
        for rec, locs in json.loads(primer.locations).iteritems():
            chr_start, chr_end = chr_ends[rec]
            new_locs += [l + chr_start for l in locs] + [chr_start, chr_end]
    return list(set(new_locs))



def binding_sites(kmer, genome_fp):
    genome = Fasta(genome_fp)
    locations = {}
    kmer = str(kmer)
    for record in genome.keys():
        seq = str(genome[record])
        locations[record] = substr_indices(kmer, seq)
        # append reversed primer locations as well
        locations[record] += substr_indices(kmer[::-1], seq)
    return locations


def _primer_bind_sites(primer, genome_fp):
    primer.locations = json.dumps(binding_sites(primer.seq, genome_fp))
    return primer


def primers_in_parallel(primers, genome_fp,
                        cores=multiprocessing.cpu_count()):
    '''
    Uses multiple processes to find the locations of all primer
    sequences in the target genome.
    '''
    updated_primers = []
    swga.message("Finding binding sites for {} primers: ".format(len(primers)))

    progressbar(0, len(primers))

    def _init_worker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    def update(_primer):
        updated_primers.append(_primer)
        progressbar(len(updated_primers), len(primers))

    pool = multiprocessing.Pool(cores, _init_worker)
    for p in primers:
        pool.apply_async(_primer_bind_sites,
                         args=(p, genome_fp),
                         callback=update)

    # Allows a keyboard interrupt to be caught whereas normally interrupts cause
    # it to hang 
    try:
        time.sleep(10)
    except KeyboardInterrupt as k:
        pool.terminate()
        pool.join()
        raise k
    else:
        pool.close()
        pool.join()
    swga.message("\n")
    return updated_primers


def chromosome_ends(genome_fp):
    '''
    Returns the locations of the starts/ends of each chromosome (record) in a
    genome where all the chromosomes are concatenated (so i.e. the 2nd genome
    start site is len(1st genome), and all indices are 0-based).
    '''
    genome = Fasta(genome_fp)
    len_so_far = 0
    chr_ends = {}
    for record in genome.keys():
        chromosome = genome[record]
        chr_len = len(chromosome)
        chr_ends[record] = [len_so_far, chr_len + len_so_far - 1]
        len_so_far += chr_len
    return chr_ends


def substr_indices(substring, string):
    '''
    Very fast way of finding overlapping substring locations in a
    (potentially large) string.
    '''
    locations = []
    start = 0
    # Assumes string is all upper-case
    substring = substring.upper()
    while True:
        start = string.find(substring, start) + 1
        if start > 0:
            locations.append(start-1)
        else:
            return locations
        
def count(substring, string):
    n = start = 0
    substring = substring.upper()
    while True:
        start = string.find(substring, start)
        if start > 0:
            n += 1
        else:
            return n

            
