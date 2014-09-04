import argparse
import ConfigParser
import os
import sys
import subprocess
import swga

def main(argv, cfg_file, quiet):
    '''Find sets of compatible primers in a primer compatibility graph.'''
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='sets',
                                   cfg_file=cfg_file)

    parser.add_argument('-i', '--input', default='-', help="""Heterodimer
    compatibility graph in DIMACS format with edges between compatible primers
    (default: stdin).""")

    parser.add_argument('-o', '--output', default='', help="""Filename to store
    results (default: stdout)""")

    parser.add_argument('-m', '--min_size', type=int, help='''Minimum size of
    primer sets. (default: %(default)s)''')

    parser.add_argument('-M', '--max_size', type=int, help='''Maximum size of
    primer sets. (default: %(default)s)''')

    parser.add_argument('-b', '--min_bg_bind_dist', type=int, help='''Minimum
    average distance between background binding sites for the primers in the
    set. (default: %(default)s)''')

    parser.add_argument('-l', '--bg_genome_len', type=int, help='''Length of
    background genome. (default: %(default)s)''')

    args = parser.parse_args(argv)
    swgahome = swga.get_swgahome()
    args.set_finder = os.path.join(swgahome, 'set_finder')

    if not os.path.isfile(args.set_finder):
        sys.stderr.write("Error: cannot find set_finder in %s.\n" % swgahome)
        exit(1)
    if not quiet:
        swga.print_cfg_file(parser.prog, cfg_file)
        swga.print_args(parser.prog, args)
    if not quiet and args.input == '-':
        swga.print_stdin_msg(parser.prog)

    find_sets(args)


def find_sets(args):
    '''
    Calls the set_finder binary with the specified options on the
    heterodimer compatibility graph and outputs valid sets for
    post-processing.

    args: Namespace object from the argument parser
    '''
    output = '> '+args.output if args.output else ''
    find_set_cmd = [args.set_finder, '-q', '-q', '-B', args.min_bg_bind_dist,
                    '-L', args.bg_genome_len, '-m', args.min_size, '-M',
                    args.max_size, '-a', '-u', '-r', 'unweighted-coloring',
                    args.input, output]
    find_set_cmd = [str(_) for _ in find_set_cmd]
    subprocess.check_call(" ".join(find_set_cmd), shell=True)

