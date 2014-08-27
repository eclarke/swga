import argparse
import ConfigParser
import os
import sys
import subprocess
import PrimerSets as ps


def main():
    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', ps.default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('find_sets'))

    parser = argparse.ArgumentParser(description="""Wrapper around set_finder to
    find sets of compatible primers.""")
    parser.set_defaults(**defaults)

    parser.add_argument('-i', '--input', default='-', help="""Heterodimer
    compatibility graph in DIMACS format with edges between compatible primers
    (default: stdin).""")

    parser.add_argument('-o', '--output', default='', help="""Filename to store
    results (default: stdout)""")

    parser.add_argument('-m', '--min_size', type=int, help='''Minimum size of
    primer sets.''')

    parser.add_argument('-M', '--max_size', type=int, help='''Maximum size of
    primer sets.''')

    parser.add_argument('-b', '--min_bg_bind_dist', type=int, help='''Minimum
    average distance between background binding sites for the primers in the
    set.''')

    parser.add_argument('-l', '--bg_genome_len', type=int, help='''Length of
    background genome.''')

    parser.add_argument('-s', '--set_finder', help='''Location of set_finder
    binary.''')

    parser.add_argument('-q', '--quiet', action='store_true', help="""Suppress
    messages (default: %(default)s)""")

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
    kwargs = vars(args)
    kwargs['output'] = '> '+kwargs['output'] if kwargs['output'] else ''
    find_set_cmd = ("{set_finder} -q -q -B {min_bg_bind_dist} -L {bg_genome_len}"
    " -m {min_size} -M {max_size} -a -u -r unweighted-coloring"
    " {input} {output}").format(**kwargs)
    subprocess.call(find_set_cmd, shell=True)

if __name__ == '__main__':
    main()
