'''
Functions used to score primer sets.
'''

from math import sqrt

def stdev(values):
    '''
    Returns the corrected sample standard deviation (with n-1 in the
    denominator rather than n) for a list of values.

    values: sequence of numeric values.
    '''
    n = float(len(values))
    mu = sum(values)/n
    return sqrt(sum((x-mu)**2 for x in values) / (n-1))


def gini(distances):
    '''
    Returns an approximation of the Gini ratio for a set of distances
    corresponding primer binding distances. The Gini ratio is between
    0 and 1, where 0 is perfectly even and 1 is very uneven, i.e. one
    primer is extremely far from the rest.

    distances: a sequence of distances. Does not need to be sorted.
    '''
    distances = sorted(distances)
    height = area = 0
    for value in distances:
        height += value
        area += height - value / 2.
    fair_area = height * len(distances) / 2
    return (fair_area - area) / fair_area


def score_seq(**kwargs):
    pass
