import argparse
import ConfigParser
import os
import sys
import PrimerSets as ps


def main(argv, cfg_file):
    defaults, _ = ps.parse_config(cfg_file, 'make_graph')
    parser = argparse.ArgumentParser(description="""Create a heterodimer
    compatibility graph from a list of primers.""",
                                     prog="swga mkgraph")
    parser.set_defaults(**defaults)

    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
    default=sys.stdin, help="""Input file where each row contains a primer, fg
    binding #, bg binding #, and fg/bg binding ratio, separated by whitespace.
    (default: stdin)""")

    parser.add_argument('-o', '--output', type=argparse.FileType('w', 0),
    default=sys.stdout, help="""Filename to store output graph (in DIMACS
    format). (default: stdout)""")

    parser.add_argument('-m', '--max_complement', type=int, help="""Max number
    of consecutive complimentary bases between two primers.""")

    parser.add_argument('-v', '--verbose', action='store_true',
    help="Display messages (default: %(default)s)")

    args = parser.parse_args(argv)
    if args.verbose and args.input.name == '<stdin>':
        sys.stderr.write("%s: Receiving input from stdin...\n" % parser.prog)
    make_graph(args)


def make_graph(args):
    '''
    Creates the heterodimer compatibility graph.

    args: Namespace object from the argument parser
    '''
    primers = ps.read_primer_file(args.input)
    arcs = ps.test_pairs(primers, args.max_hetdimer_bind)
    ps.write_graph(primers, arcs, args.output)
