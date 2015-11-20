# -*- coding: utf-8 -*-
"""score.py

Functions for retrieving and scoring primer sets.

"""

import importlib

import locate
import stats

from swga import warn


def read_set_finder_line(line):
    '''
    Reads a line in the format [primer_id1,primer_id2,... weight] and returns
    a tuple (size, weight, [primer_id1, primer_id2, ...]).
    '''
    primer_set, weight = line.strip('\n').split(' ')
    primer_set = [int(_) for _ in primer_set.split(',')]
    return (primer_set, float(weight))


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


def default_score_set(expression, primer_set, primer_locs, max_dist, bg_dist_mean):
    # Calculate various metrics
    binding_distances = stats.seq_diff(primer_locs)
    namespace = {
        'set_size': len(primer_set),
        'fg_dist_mean': stats.mean(binding_distances),
        'fg_dist_std': stats.stdev(binding_distances),
        'fg_dist_gini': stats.gini(binding_distances),
        'bg_dist_mean': bg_dist_mean,
        'fg_max_dist': max_dist,
        '__builtins__': None}
    permitted_var_str = ", ".join(
        [key for key in namespace.keys() if key is not "__builtins__"])
    score = None
    try:
        score = eval(expression, namespace, {'__builtins__': {}})
    except NameError as e:
        raise NameError(
            e.message +
            '. Permitted variables are %s. Refer to README or docs for help.'
            % permitted_var_str)
    del namespace['__builtins__']
#    print_primer_set(primer_set, [score, namespace], output_handle)
    return score, namespace


def calculate_bg_dist_mean(primers, bg_length):
    total_bg_freq = sum(p.bg_freq for p in primers)
    if total_bg_freq == 0:
        warn(
            "No primers appear in the background genome: "
            "bg_dist_mean set as infinite")
        bg_dist_mean = float('Inf')
    else:
        bg_dist_mean = float(
            bg_length) / sum(p.bg_freq for p in primers)
    return bg_dist_mean


def score_set(primers, max_fg_bind_dist, bg_dist_mean,
              chr_ends, score_fun, interactive=False):

    binding_locations = locate.linearize_binding_sites(
        primers, chr_ends)
    max_dist = max(stats.seq_diff(binding_locations))

    # If it's not a user-supplied set and it's not passing the filter,
    # abort immediately
    if not interactive and max_dist > max_fg_bind_dist:
        return False, {}, max_dist

    set_score, variables = score_fun(
        primer_set=primers,
        primer_locs=binding_locations,
        max_dist=max_dist,
        bg_dist_mean=bg_dist_mean)

    return set_score, variables, max_dist
