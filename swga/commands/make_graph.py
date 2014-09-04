import argparse
import ConfigParser
import os
import sys
import swga


def main(argv, cfg_file, quiet):
    '''Create a heterodimer compatibility graph from a list of primers.'''
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='mkgraph',
                                   cfg_file=cfg_file)

    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
    default=sys.stdin, help="""Input file where each row contains a primer, fg
    binding #, bg binding #, and fg/bg binding ratio, separated by whitespace.
    (default: stdin)""")

    parser.add_argument('-o', '--output', type=argparse.FileType('w', 0),
    default=sys.stdout, help="""Filename to store output graph (in DIMACS
    format). (default: stdout)""")

    parser.add_argument('-m', '--max_complement', type=int, help="""Max number
    of consecutive complimentary bases between two primers.""")

    args = parser.parse_args(argv)
    if not quiet and args.input.name == '<stdin>':
        swga.print_stdin_msg(parser.prog)

    make_graph(args)


def make_graph(args):
    '''
    Creates the heterodimer compatibility graph.

    args: Namespace object from the argument parser
    '''
    primers = swga.read_primer_file(args.input)
    arcs = swga.test_pairs(primers, args.max_hetdimer_bind)
    swga.write_graph(primers, arcs, args.output)
