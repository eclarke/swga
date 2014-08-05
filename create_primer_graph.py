import PrimerSets as ps
from argparse import ArgumentParser
from os.path import isfile
import sys

def main():
    
    opts = ps.read_config_file(ps.default_config_file)
    if not isfile(ps.default_config_file):
        opts.set('primer_filters', 'max_htdmr_bind', '0')
        sys.sterr.write(ps.opts_errstr)
    
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

    primers = ps.read_primers(args['primer_file'])
    arcs = ps.test_pairs(primers, args['max_htdmr_bind'])
    ps.write_graph(primers, arcs, args['output'])


if __name__ == '__main__':
        main()
