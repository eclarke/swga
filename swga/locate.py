# -*- coding: utf-8 -*-
"""locate.py

This module contains functions for locating primers in genomes (i.e. finding
small substrings in larger strings).

"""
import json
from pyfaidx import Fasta


def revcomp(s):
    r = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G'
    }
    rc = [r[i] for i in s][::-1]
    return "".join(rc)


def linearize_binding_sites(primers, chr_ends):
    '''
    Modifies the primer binding site locations as if they were positions on a
    linear genome composed of all the chromosomes concatenated together. This
    allows us to compute distances between primer binding sites correctly.
    '''
    new_locs = []
    for primer in primers:
        for rec, locs in primer.locations.iteritems():
            chr_start, chr_end = chr_ends[rec]
            new_locs += [l + chr_start for l in locs] + [chr_start, chr_end]
    new_locs = list(set(new_locs))
    if new_locs == []:
        raise ValueError("Binding sites for primers not found!")
    return list(set(new_locs))


def binding_sites(kmer, genome_fp):
    genome = Fasta(genome_fp)
    locations = {}
    kmer = str(kmer)
    for record in genome.keys():
        seq = str(genome[record])
        locations[record] = substr_indices(kmer, seq)
        # append reversed primer locations as well
        locations[record] += substr_indices(revcomp(kmer), seq)
    if locations == {}:
        raise ValueError(
            "No locations for {} found in fg genome!".format(kmer))
    return locations


def _primer_bind_sites(primer, genome_fp):
    primer.locations = json.dumps(binding_sites(primer.seq, genome_fp))
    return primer


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
