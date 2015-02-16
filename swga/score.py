# -*- coding: utf-8 -*-
"""score.py

Functions for retrieving and scoring primer sets.

"""
import stats
import importlib
import json

def read_set_finder_line(line):
    '''
    Reads a line in the format [primer_id1,primer_id2,... weight] and returns
    a tuple (size, weight, [primer_id1, primer_id2, ...]).
    '''
    primer_set, weight = line.strip('\n').split(' ')
    primer_set = [int(_) for _ in primer_set.split(',')]
    return (primer_set, float(weight))


def linearize_binding_sites(primers, chr_ends):
    '''
    Modifies the primer binding site locations as if they were positions on a
    linear genome composed of all the chromosomes concatenated together. This
    allows us to compute distances between primer binding sites correctly.
    '''
    new_locs = []
    for primer in primers:
        for rec, locs in json.loads(primer.locations).iteritems():
            print rec
            assert rec in chr_ends.keys()
            chr_start, chr_end = chr_ends[rec]
            new_locs += [l + chr_start for l in locs] + [chr_start, chr_end]
    return new_locs


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
    string would be swga:score:get_user_fun.
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
#    print_primer_set(primer_set, [score, namespace], output_handle)
    return score, namespace


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


