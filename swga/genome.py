import mmap
import cPickle
import gzip
import signal
import time
import sys
import multiprocessing
from contextlib import closing

from .core import progressbar

def check_if_flattened(fasta_fp):
    '''Checks if file only has one line'''
    with open(fasta_fp) as fasta:
        for i, _ in enumerate(fasta):
            if i > 0:
                return False
        return True


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
        sys.stderr.write("Finding primer locations...\n")
        progressbar(0, len(primers))

    def update_locations(loc):
        primer = loc[0]
        p_locs = loc[1]
        locations[primer.id] = {'primer': primer,
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


def save_locations(locations, filename, verbose=False):
    '''Saves primer binding locations to a gzipped pickle file.'''
    with gzip.GzipFile(filename, 'w') as dest:
        cPickle.dump(locations, dest)
        if verbose:
            sys.stderr.write("Locations stored in %s\n" % dest.name)


def load_locations(filename):
    '''Loads primer binding locations from gzipped pickled file.'''
    try:
        with gzip.GzipFile(filename, 'r') as f:
            try:
                return cPickle.load(f)
            except (IOError, cPickle.UnpicklingError):
                raise IOError("Cannot read primer locations from file '%s'" % filename)
    except TypeError:
        raise ValueError("No primer location file specified.")
