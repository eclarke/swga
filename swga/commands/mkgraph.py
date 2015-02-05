import swga
import swga.primers
import swga.graph as graph
from swga.commands import Command


def main(argv, cfg_file):
    cmd = Command('mkgraph', cfg_file = cfg_file)
    cmd.parse_args(argv)
    make_graph(**cmd.args)


def make_graph(input,
               output,
               max_hetdimer_bind):
    '''
    Creates a heterodimer compatibility graph.
    '''
    primers = swga.primers.read_primer_file(input)
    arcs = graph.test_pairs(primers, max_hetdimer_bind)
    graph.write_graph(primers, arcs, output)
