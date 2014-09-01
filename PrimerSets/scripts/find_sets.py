import argparse
import ConfigParser
import os
import sys
import subprocess
import PrimerSets as ps


def main(argv, cfg_file):
    config = ConfigParser.SafeConfigParser()
    if not cfg_file:
        cfg_file = os.environ.get('swga_params', ps.default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('find_sets'))

    parser = argparse.ArgumentParser(description="""Wrapper around set_finder to
    find sets of compatible primers.""", prog="swga sets")
    parser.set_defaults(**defaults)

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

    parser.add_argument('-s', '--set_finder', help='''Location of set_finder
    binary. (default: %(default)s)''')

    parser.add_argument('-q', '--quiet', action='store_true', help="""Suppress
    messages (default: %(default)s)""")

    args = None
    if argv:
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()
    if not args.quiet and args.input == '-':
        sys.stderr.write("Receiving input from stdin...\n")
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

if __name__ == '__main__':
    main(None, None)
