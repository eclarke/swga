import itertools
from .primers import Primer

def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of edges if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    edges = []
    for p1, p2 in itertools.combinations(starting_primers, 2):
        if (p1.seq not in p2.seq) and (p2.seq not in p1.seq):
            if max_consecutive_binding(p1.seq, p2.seq) <= max_binding:
                edges.append([p1._id, p2._id])
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
            weight = primer.bg_freq if primer.bg_freq > 0 else 1
            file_handle.write('n {} {}\n'.format(primer._id,
                                                 weight))

        except AttributeError:
            raise ValueError("Primers must be of the form {}".format(type(Primer)))
    for edge in edges:
        try:
            file_handle.write('e {} {}\n'.format(edge[0], edge[1]))
        except IndexError:
            raise ValueError("Edges must be specified as a list with"+
            "two elements. Invalid edge: {}".format(edge))
