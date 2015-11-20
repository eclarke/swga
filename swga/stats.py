"""stats.py

Functions used to score primer sets.

"""

from math import sqrt


def mean(values):
    """Calculate the arithmetic mean."""
    return float(sum(values)) / len(values) if len(values) > 0 else float('nan')


def stdev(values):
    """Calculate the corrected sample standard deviation."""
    n = float(len(values))
    mu = sum(values) / n
    return sqrt(sum((x - mu)**2 for x in values) / (n - 1))


def seq_diff(seq):
    """Calculate the sequential differences of a set of numbers.

    :param seq: a set of numbers, will be sorted before sequential differences
    """
    seq.sort()
    diffs = []
    for i in xrange(len(seq) - 1):
        diff = seq[i + 1] - seq[i]
        assert diff >= 0
        diffs.append(diff)
    return diffs


def gini(distances):
    """Calculate an approximation of the Gini coefficient for a set of values.

    The Gini coefficient represents the evenness of a set of values, and is in
    the range [0,1], where 0 is perfectly even and 1 is perfectly uneven.

    In practice, this means that sets with low Gini coefficients have more even
    binding frequencies and are more likely to not favor highly-repeated regions
    or other aberrations in the target genome.

    :param distances: a vector of primer binding distances (from each other)
    """
    distances.sort()
    height = area = 0
    for value in distances:
        height += value
        area += height - value / 2.
    fair_area = height * len(distances) / 2
    return (fair_area - area) / fair_area


def lorenz(distances):
    """Calculate the Lorenz curve for a set of values.

    The Lorenz curve is the cumulative total value of a set of values ordered
    by smallest to largest. The more this deflects from the line of perfect
    equality, the higher the relative inequality of the values (the Gini
    coefficient. These numbers are provided to plot the Lorenz curve for
    visualization purposes.

    For more information, see https://en.wikipedia.org/wiki/Lorenz_curve.

    :param distances: a vector of primer binding distances (from each other)
    """
    distances.sort()
    height = 0
    lz = []
    for dist in distances:
        height += dist
        lz.append(height)

    # Normalize output
    _max = float(max(lz))
    lz = [l / _max for l in lz]

    return lz
