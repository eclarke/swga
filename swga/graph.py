import itertools
from swga.kmers import max_sequential_nt


def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of edges if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    edges = []
    for p1, p2 in itertools.combinations(starting_primers, 2):
        if (p1.seq not in p2.seq) and (p2.seq not in p1.seq):
            if max_sequential_nt(p1.seq, p2.seq) <= max_binding:
                edges.append([p1._id, p2._id])
    return edges


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
        weight = primer.bg_freq if primer.bg_freq > 0 else 1
        file_handle.write('n {} {}\n'.format(primer._id, weight))
    for edge in edges:
        try:
            file_handle.write('e {} {}\n'.format(edge[0], edge[1]))
        except IndexError:
            raise ValueError("Edges must be specified as a list with" +
                             "two elements. Invalid edge: {}".format(edge))
