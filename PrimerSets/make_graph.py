import argparse
import ConfigParser
import os
import sys
import PrimerSets as ps


def main():
    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', ps.default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('make_graph'))

    parser = argparse.ArgumentParser(description="""Create a heterodimer
    compatibility graph from a list of primers.""")
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

    parser.add_argument('-q', '--quiet', action='store_true',
    help="Suppress messages (default: %(default)s)")

    args = parser.parse_args()
    if not args.quiet and args.input.name == '<stdin>':
        sys.stderr.write("Receiving input from stdin...\n")
    make_graph(args)


def make_graph(args):
    '''
    Creates the heterodimer compatibility graph.

    args: Namespace object from the argument parser
    '''
    primers = ps.read_primer_file(args.input)
    arcs = ps.test_pairs(primers, args.max_hetdimer_bind)
    ps.write_graph(primers, arcs, args.output)


if __name__ == '__main__':
    main()
