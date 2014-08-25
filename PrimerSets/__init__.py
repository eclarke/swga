import re
import sys
import time
import mmap
import signal
import itertools
import multiprocessing
from collections import namedtuple
from contextlib import closing

# Definitions
Primer = namedtuple('Primer', 'id, seq, bg_freq, fg_freq, ratio')
default_config_file = 'parameters.cfg'

# Functions
def parse_primer(string, line_no=0):
    '''
    Takes a line from a tab- or space-delimited file where each row specifies
    a primer sequence, fg binding count, bg binding count, and fg/bg binding
    count ratio.

    Returns: A Primer object (or raises ValueError if it cannot parse the line)
    '''
    seq, fg_freq, bg_freq, ratio = re.split(r'[ \t]+', string.strip('\n'))
    return Primer(line_no, seq, int(bg_freq), int(fg_freq), float(ratio))


def read_primer_file(file_handle, echo_input=False, verbose=False):
    '''
    Calls parse_primer() on each line of the input file. If a malformed line is
    found, will skip parsing and (optionally) output a warning message.

    Arguments:
    file_handle: an open file handle to read from
    echo_input:  if True, echo each line of the file to stdout
    verbose:     if True, warn when parsing a line fails

    Returns: A list of Primer objects
    '''
    primers = []
    for i, line in enumerate(file_handle):
        try:
            primers.append(parse_primer(line, i))
            if echo_input:
                sys.stdout.write(line)
        except ValueError as e:
            if verbose:
                sys.stderr.write("Cannot parse line %i (reason: %s), "\
                "skipping...\n" % (i, e.message))
            continue
    return primers

## Heterodimer graph functions

def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of edges if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    edges = []
    for p1, p2 in itertools.combinations(starting_primers, 2):
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


## Genome binding location functions

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


def mp_find_primer_locations(primers, genome_fp, cores, verbose):
    '''
    Uses multiple processes to find the locations of all primer
    sequences in the target genome.
    '''
    locations = {}
    if verbose:
        progressbar(0, len(primers))
    def update_locations(loc):
        primer = loc[0]
        p_locs = loc[1]
        locations[primer.id] = {'seq': primer.seq,
                                'loc': p_locs}
        if verbose:
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
    if verbose:
        sys.stderr.write('\n')
    return locations


def progressbar(i, length):
    if i >= 1:
        i = i/(length*1.0)
    sys.stderr.write('\r[%-20s] %-3d%%' % ('='*int(round(i*20)), i*100))
    sys.stderr.flush()
