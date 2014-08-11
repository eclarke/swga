from math import sqrt

def stdev(seq):
    n = float(len(seq))
    mu = sum(seq)/n
    return sqrt(sum((x-mu)**2 for x in seq) / (n-1))

def gini(seq):
    sorted_seq = sorted(seq)
    height = area = 0
    for value in sorted_seq:
        height += value
        area += height - value / 2.
    fair_area = height * len(sorted_s) / 2
    return (fair_area - area) / fair_area
