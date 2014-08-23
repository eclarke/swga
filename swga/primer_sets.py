# primer_sets.py
# Functions to read in a list of primers, test primer pairs for
# heterodimer-ness, and write out a DIMACS graph file for use
# downstream in cliquer.
#
# Erik Clarke - ecl@mail.med.upenn.edu

from itertools import combinations
from argparse import ArgumentParser
from collections import namedtuple
from ConfigParser import SafeConfigParser
import mmap
import multiprocessing
import sys
import os
import re
import itertools
import signal
import time
from contextlib import closing

Primer = namedtuple('Primer', 'id, seq, bg_freq, fg_freq')

default_config_file = 'parameters.cfg'

# formatted by downstream modules
opts_errstr = """
--WARNING: Cannot find default config file, specified as '{}'.--
Ensure all values are specified on the command line, or set and export
the swga_params environment variable. Unspecified parameters set to
None or 0. \n
"""

def read_config_file(filename):
    parser = SafeConfigParser()
    parser.read(filename)
    if not parser.has_section('primer_filters'):
        parser.add_section('primer_filters')
    if not parser.has_section('set_opts'):
        parser.add_section('set_opts')
    if not parser.has_section('locations'):
        parser.add_section('locations')
    return parser


def write_graph(primers, edges, file_handle):
    '''
    Writes graph in DIMACS graph format, specified in the cliquer user manual.
    See http://users.tkk.fi/~pat/cliquer.html
    
    An edge is a list of the form [first_node, second_node]
    "edges" is a list of edges
    '''

    if type(file_handle) is not file:
        raise ValueError(("file_handle must be a file, not "
                          "a {}").format(type(file_handle))) 

    file_handle.write('p sp {} {}\n'.format(len(primers), len(edges)))
    for primer in primers:
        try:
            file_handle.write('n {} {}\n'.format(primer.id,
                                                 primer.bg_freq))
        except AttributeError:
            raise ValueError("Primers must be of the form {}".format(type(Primer)))
    for edge in edges:
        try:
            file_handle.write('e {} {}\n'.format(edge[0], edge[1]))
        except IndexError:
            raise ValueError("Edges must be specified as a list with"+
            "two elements. Invalid edge: {}".format(edge)) 


def read_primers(file_handle):
    '''
    Reads in a space-delimited file where the first column is the
    primer sequence, second column is the foreground genome bind
    count, and third is the background genome binding count. 
    Returns a list of Primer objects.
    '''
    try:
        primers = []
        for i, line in enumerate(file_handle):
            seq, fg_freq, bg_freq, ratio = line.strip('\n').split(' ')
            primers.append(Primer(len(primers)+1, seq,
                                  int(bg_freq), int(fg_freq)))
    except ValueError as err:
        
        raise err
    return primers


def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of edges if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    edges = []
    for p1, p2 in combinations(starting_primers, 2):
        if max_consecutive_binding(p1.seq, p2.seq) <= max_binding:
            edges.append([p1.id, p2.id])
    return edges


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


def find_locations(substring, string):
    '''
    Very fast way of finding overlapping substring locations in a
    (potentially large) string.
    '''
    locations = []
    start = 0
    while True:
        start = string.find(substring, start) + 1
        if start > 0:
            locations.append(start-1)
        else:
            return locations


def find_primer_locations(primer, genome_fp):
    '''
    Returns the matches in the genome of the primer, reverse primer,
    and the chromosome end points. May be imprecise as the chromosome
    start marker is occupied by a char (>), so possible off-by-one errors.
    '''
    with open(genome_fp) as f:
        with closing(mmap.mmap(f.fileno(), 0,
                               access=mmap.ACCESS_COPY)) as genome:
            record_locations = find_locations('>', genome)
            if len(record_locations) == 0:
                raise ValueError('No records found in genome!')
            primer_locations = find_locations(primer.seq, genome)
            rev_primer_locations = find_locations(primer.seq[::-1], genome)
            return (primer, sorted(primer_locations + rev_primer_locations + record_locations))


def _init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def mp_find_primer_locations(primers, genome_fp,
                             cores=multiprocessing.cpu_count(),
                             chatty=True):
    '''
    Uses multiple processes to find the locations of all primer
    sequences in the target genome. 
    '''
    locations = {}
    if chatty:
        progressbar(0, len(primers))
    def update_locations(loc):
        primer = loc[0]
        p_locs = loc[1]
        locations[primer.id] = {'seq': primer.seq,
                                'loc': p_locs}
        if chatty:
            progressbar(len(locations.keys()), len(primers))

    
    pool = multiprocessing.Pool(cores, _init_worker)
    for primer in primers:
        pool.apply_async(find_primer_locations,
                         args=(primer, genome_fp),
                         callback=update_locations)

    # it's unclear to me why this works, but it allows a keyboard
    # interrupt to be caught whereas normally interrupts cause it
    # to hang
    try:
        time.sleep(10)
    except KeyboardInterrupt as k:
        pool.terminate()
        pool.join()
        raise k
    else:
        pool.close()
        pool.join()
    if chatty:
        sys.stderr.write('\n')
    return locations


def fg_bind_distances(setline, primer_locations, stat_func):
    pset_line = setline.strip('\n').split(' ')
    psize = pset_line[0]
    pweight = pset_line[1]
    primer_set = [int(_) for _ in pset_line[2::]]
    locations = sum([primer_locations[primer]['loc'] for primer in
                    primer_set], [])
    primers = [primer_locations[primer]['seq'] for primer in \
               primer_set] 
    stdev = stat_func(locations)
    max_dist = max_seq_diff(sorted(locations))
    return (primer_set, primers, max_dist, stdev)


def max_seq_diff(seq):
    '''
    Calculates the sequential difference along a sorted sequence of
    integers and returns the max value.
    '''
    max_diff = 0
    for i in range(len(seq)-1):
        diff = seq[i+1] - seq[i]
        assert diff >= 0
        if diff > max_diff:
            max_diff = diff
    return max_diff


def stdev(seq):
    from math import sqrt
    n = float(len(seq))
    mu = sum(seq)/n
    return sqrt(sum((x-mu)**2 for x in seq) / (n-1))
    

def progressbar(i, length):
    if i >= 1:
        i = i/(length*1.0)
    sys.stderr.write('\r[%-20s] %-3d%%' % ('='*int(round(i*20)), i*100))
    sys.stderr.flush()

    
class DefaultValueError(Exception):
    def __init__(self, missing_arg):
        super(DefaultValueError, self).__init__(missing_arg+" not specified and no"+
                                            " default found in config file.")

class MaxSetsSignal(Exception):
    pass
