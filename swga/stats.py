"""stats.py

Functions used to score primer sets.

"""

from math import sqrt


def mean(values):
    '''Returns the arithmetic mean.'''
    return float(sum(values))/len(values) if len(values) > 0 else float('nan')


def stdev(values):
    '''
    Returns the corrected sample standard deviation (with n-1 in the
    denominator rather than n) for a list of values.

    values: sequence of numeric values.
    '''
    n = float(len(values))
    mu = sum(values)/n
    return sqrt(sum((x-mu)**2 for x in values) / (n-1))


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


def gini(distances):
    '''
    Returns an approximation of the Gini ratio for a set of distances
    corresponding primer binding distances. The Gini ratio is between
    0 and 1, where 0 is perfectly even and 1 is very uneven, i.e. one
    primer is extremely far from the rest.

    distances: a sequence of distances. Does not need to be sorted.
    '''

    distances.sort()
    height = area = 0
    for value in distances:
        height += value
        area += height - value / 2.
    fair_area = height * len(distances) / 2
    return (fair_area - area) / fair_area


def lorenz(distances):
    '''Returns the Lorenz curve of the inter-primer distances'''

    distances.sort()
    height = 0
    lz = []
    for dist in distances:
        height += dist
        lz.append(height)

    # Normalize output
    _max = float(max(lz))
    lz = [l/_max for l in lz]

    return lz


