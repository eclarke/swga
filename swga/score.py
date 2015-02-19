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
    string would be swga.score:get_user_fun.
    '''
    try:
        module, fun = spec_str.split(':')
    except ValueError:
        raise ValueError("Invalid function specification string. Must have the "
        "format modulename.possible_submodule:function_name""")
    module = importlib.import_module(module)
    return getattr(module, fun)


def default_score_set(expression, primer_set, primer_locs, max_dist, bg_ratio):
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


