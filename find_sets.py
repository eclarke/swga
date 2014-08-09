#!/usr/bin/env python

import primer_sets as ps
from argparse import ArgumentParser
import os
import sys
import subprocess

def main():
    config_file = os.environ.get('swga_params', ps.default_config_file)
    opts = ps.read_config_file(config_file)
    if not os.path.isfile(ps.default_config_file):
        opts.set('set_opts','min_size', '0')
        opts.set('set_opts','max_size', '0')
        opts.set('set_opts','min_bg_bind_dist', '0')
        opts.set('set_opts','bg_genome_len', '0')
        opts.set('locations','set_finder', '0')
        sys.stderr.write(ps.opts_errstr.format(config_file))
    
    description_string = '''Reads in a primer graph and finds
    compatible primer sets of the specified size range and background
    binding frequency. Writes results to specified output file or
    stdout if output unspecified.'''

    parser = ArgumentParser(description=description_string)

    parser.add_argument('-m', '--min_size', type=int,
                        default=opts.getint('set_opts', 'min_size'),
                        help='''Minimum size of primer sets.''')
    parser.add_argument('-M', '--max_size', type=int,
                        default=opts.getint('set_opts', 'max_size'),
                        help='''Maximum size of primer sets.''')
    parser.add_argument('-b', '--min_bg_bind_dist', type=int,
                        default=opts.getint('set_opts',
                                            'min_bg_bind_dist'),
                        help='''Minimum average distance between
                        background binding sites for the primers in
                        the set.''')
    parser.add_argument('-l', '--bg_genome_len', type=int,
                        default=opts.getint('set_opts',
                                             'bg_genome_len'),
                        help='''Length of background genome.''')

    parser.add_argument('-s', '--set_finder',
                        default=opts.get('locations',
                                          'set_finder'),
                        help='''Location of set_finder binary.''')

    parser.add_argument('-g', '--graph_file', action='store',
                        help='''Primer graph in DIMACS format where
                        edges are between compatible primers. Default
                        to stdin if unspecified.''',
                        default=None)

    parser.add_argument('-o', '--output', action='store',
                        help='''Where to store results. Default to
                        stdout if unspecified.''', default=None)

    args = vars(parser.parse_args())

    command = """{set_finder} -q -q -B {min_bg_bind_dist} -L {bg_genome_len} \
    -m {min_size} -M {max_size} -a -u -r unweighted-coloring \
    {graph_file} {output}""".format(
    set_finder=args['set_finder'],
    min_bg_bind_dist=args['min_bg_bind_dist'],
    bg_genome_len=args['bg_genome_len'],
    min_size=args['min_size'],
    max_size=args['max_size'],
    graph_file=args['graph_file'] if args['graph_file'] else '-',
    output='> ' + args['output'] if args['output'] else '')
    subprocess.call(command, shell=True)


if __name__ == "__main__":
    main()
    
