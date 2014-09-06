import re
import os
import sys
import gzip
import time
import mmap
import stats
import errno
import signal
import cPickle
import argparse
import itertools
import importlib
import ConfigParser
import multiprocessing

from collections import namedtuple
from contextlib import closing

# Definitions
Primer = namedtuple('Primer', 'id, seq, bg_freq, fg_freq, ratio')
default_config_file = 'parameters.cfg'

# Functions
def get_swgahome():
    swgahome = os.environ.get('SWGAHOME')
    if not swgahome:
        raise ValueError("SWGAHOME not set, cannot find home directory for SWGA scripts.")
    if not os.path.isabs(swgahome):
        raise ValueError("SWGAHOME cannot be a relative path. Make SWGAHOME an absolute path and try again.")
    return swgahome

def mkdirp(path):
    '''Simulates 'mkdir -p': creates a directory unless it already exists'''
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def basic_cmd_parser(description, cmd_name, cfg_file):
    defaults, _ = parse_config(cfg_file, cmd_name)
    print defaults, _
    parser = argparse.ArgumentParser(description=description, prog='swga '+cmd_name)
    parser.set_defaults(**defaults)
    return parser

def print_stdin_msg(prog_name):
    sys.stderr.write("{}: receiving input from stdin...\n".format(prog_name))

def print_args(prog_name, args):
    # argstr = json.dumps(args, sort_keys=True, indent=2, separators=(',', ': '))
    sys.stderr.write("{} parameters: {}\n".format(prog_name, str(vars(args))))

def print_cfg_file(prog_name, cfg_file):
    sys.stderr.write("{} config file: {}\n".format(prog_name, os.path.abspath(cfg_file)))

def parse_config(cfg_file, section):
    '''
    Parses a config file and returns a dictionary of the values found
    in the specified section, along with the ConfigParser itself
    '''
    config = ConfigParser.SafeConfigParser()
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items(section))
    return defaults, config


def parse_primer(string, line_no=1):
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
            primers.append(parse_primer(line, i+1))
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

## Set processing functions
def read_set_finder_line(line):
    '''
    Reads a line in the format [size weight primer_id1,primer_id2,...] and returns
    a tuple (size, weight, [primer_id1, primer_id2, ...]).
    '''
    primer_set, weight = line.strip('\n').split(' ')
    primer_set = [int(_)-1 for _ in primer_set.split(',')]
    return (primer_set, float(weight))


def get_primers_from_ids(primer_ids, primer_store):
    '''
    Retrieves the Primer object for each id in a list from the stored locations.

    Arguments:
    primer_ids: a list of primer ids (integers)
    primer_store: A dict of the form {primer_id: {'primer':Primer, 'loc':[locations]}}

    Returns: a list of Primers
    '''
    return [primer_store[primer]['primer'] for primer in primer_ids]


def get_primer_locations(primer_ids, primer_store):
    '''
    Retrieves the primer binding locations for each id in a list from the stored
    binding locations.

    Arguments:
    primer_ids: a list of primer ids (integers)
    primer_store: A dict of the form {primer_id: {'primer':Primer, 'loc':[locations]}}

    Returns: a list with all the binding sites of the primers in a set, aggregated
    '''
    # Aggregates all the locations into one list
    return sum([primer_store[primer]['loc'] for primer in primer_ids], [])


def seq_diff(seq):
    '''
    Returns the sequential differences along a sorted sequence of numbers.
    If the sequence is not already sorted, it will sort it first.
    '''
    seq.sort()
    diffs = []
    for i in xrange(len(seq)-1):
        diff = seq[i+1] - seq[i]
        assert diff >= 0
        diffs.append(diff)
    return diffs


def get_user_fun(spec_str):
    '''
    Parses a string to get a function from a module. The string format is simply
    modulename.possible_submodule:function_name. For instance, this function's
    string would be PrimerSets:get_user_fun.
    '''
    try:
        module, fun = spec_str.split(':')
    except ValueError:
        raise ValueError("Invalid function specification string. Must have the "
        "format modulename.possible_submodule:function_name""")
    module = importlib.import_module(module)
    return getattr(module, fun)


def default_score_set(expression, primer_set, primer_locs, max_dist, bg_ratio,
    output_handle):
    # Calculate various metrics
    binding_distances = seq_diff(primer_locs)
    namespace = {
        'set_size': len(primer_set),
        'fg_dist_mean': stats.mean(binding_distances),
        'fg_dist_std': stats.stdev(binding_distances),
        'fg_dist_gini': stats.gini(binding_distances),
        'bg_ratio': bg_ratio,
        'fg_max_dist': max_dist,
        '__builtins__': None}
    permitted_var_str = ", ".join([key for key in namespace.keys() if key is not "__builtins__"])
    score = None
    try:
        score = eval(expression, namespace, {'__builtins__': {}})
    except NameError as e:
        raise NameError(e.message + '. Permitted variables are %s. Refer to README or docs for help.' % permitted_var_str)
    del namespace['__builtins__']
    print_primer_set(primer_set, [score, namespace], output_handle)


def print_primer_set(primers, other_vals, output_handle):
    '''
    Writes the primer sequences (joined by commas) and the other_vals (separated
    by tabs) to the specified output.
    Example output:
    AAATTT,GGGCCC,ATGCATGC  val1    val2    val3
    '''
    primers_str = ",".join([primer.seq for primer in primers])
    line = "\t".join([primers_str] + [str(_) for _ in other_vals])
    output_handle.write(line+'\n')


def progressbar(i, length):
    if i >= 1:
        i = i/(length*1.0)
    sys.stderr.write('\r[%-20s] %-3d%%' % ('='*int(round(i*20)), i*100))
    sys.stderr.flush()
