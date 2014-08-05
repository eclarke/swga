# primer_cliques.py
# Functions to read in a list of primers, test primer pairs for
# heterodimer-ness, and write out a DIMACS graph file for use
# downstream in cliquer.
#
# Erik Clarke - ecl@mail.med.upenn.edu

from itertools import combinations
from argparse import ArgumentParser
from collections import namedtuple
import sys

Primer = namedtuple('Primer', 'id, seq, bg_freq, fg_freq')

def write_graph(primers, arcs, fname):
    '''
    Writes graph in DIMACS graph format, specified in the cliquer user manual.
    See http://users.tkk.fi/~pat/cliquer.html
    
    An arc is simply a list of the form [first_node, second_node]
    "arcs" is a list of arcs
    '''
    num_nodes = len(primers)
    with open(fname, 'w') as output:
        output.write('p sp {} {}\n'.format(num_nodes, len(arcs)))
        for primer in primers:
            output.write('n {} {}\n'.format(primer.id, primer.bg_freq))
        for arc in arcs:
            output.write('e {} {} \n'.format(arc[0], arc[1]))


def read_primers(primer_fp):
    '''
    Reads in a tab-delimited file where the first column is the primer
    sequence, second column is the foreground genome bind count, and
    third is the background genome binding count.
    Returns a list of Primer objects.
    '''
    with open(primer_fp) as primer_file:
        primers = []
        try:
            for i, line in enumerate(primer_file):
                seq, fg_freq, bg_freq, ratio = line.strip('\n').split(' ')
                primers.append(Primer(len(primers)+1, seq,
                                      int(bg_freq), int(fg_freq)))
        except ValueError as err:
            sys.stderr.write("Invalid primer file format.\n")
            raise err
    return primers


def write_primers(primers, fname):
    with open(fname, 'w') as output:
        output.write('ID\tSEQ\tBG_FREQ\tFG_FREQ\n')
        [output.write("\t".join([str(p.id), p.seq, str(p.bg_freq),
                                str(p.fg_freq)])+'\n') for p in primers]


def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of "arcs" if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    arcs = []
    for p1, p2 in combinations(starting_primers, 2):
        if max_consecutive_binding(p1.seq, p2.seq) < max_binding:
            arcs.append([p1.id, p2.id])
    return arcs


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
