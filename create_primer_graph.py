#!/usr/bin/env python

import primer_sets as ps
from argparse import ArgumentParser
import os
import sys

def main():
    config_file = os.environ.get('swga_params', ps.default_config_file)
    opts = ps.read_config_file(config_file)
    if not os.path.isfile(ps.default_config_file):
        opts.set('primer_filters', 'max_htdmr_bind', '0')
        sys.stderr.write(ps.opts_errstr.format(config_file))
    
    description_string = '''Reads in primers from a tab-delimited
    file, performs heterodimer checks and then writes results in
    DIMACS format for use with cliquer.c.''' 

    parser = ArgumentParser(description = description_string)

    parser.add_argument('-m', '--max_htdmr_bind', help='''Max number of
    consecutive complimentary bases allowed in the heterodimer
    filter''', type=int, default=opts.getint('primer_filters',
                                             'max_htdmr_bind')) 
    
    parser.add_argument('-o', '--output', help='''Filename to store the
    DIMACS output graph. If unspecified, writes to stdout.''', default=sys.stdout)

    parser.add_argument('-p', '--primer_file', action='store', help='''List of
    primers (first col), number of foreground binding sites (second
    col), number of background binding sites (third col), and fg/bg
    ratio (fourth col). If unspecified, reads from stdin.''', default=sys.stdin)

    args = vars(parser.parse_args())
    primers = []
    if type(args['primer_file']) == str:
        with open(args['primer_file']) as infile:
            primers = ps.read_primers(infile)
    else:
        primers = ps.read_primers(args['primer_file'])

    arcs = ps.test_pairs(primers, args['max_htdmr_bind'])

    if type(args['output']) == str:
        with open(args['output'], 'w') as outfile:
            ps.write_graph(primers, arcs, outfile)
    else:
        ps.write_graph(primers, arcs, args['output'])

if __name__ == '__main__':
        main()
