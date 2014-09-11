import stats
import importlib


def read_set_finder_line(line):
    '''
    Reads a line in the format [size weight primer_id1,primer_id2,...] and returns
    a tuple (size, weight, [primer_id1, primer_id2, ...]).
    '''
    primer_set, weight = line.strip('\n').split(' ')
    primer_set = [int(_) for _ in primer_set.split(',')]
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


