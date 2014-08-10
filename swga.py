#!/usr/bin/env python

import primer_sets as ps
import shlex
import sys
import argparse
import ConfigParser
import subprocess
import cPickle
import gzip
import multiprocessing
from signal import signal, SIGPIPE, SIG_DFL, SIGTERM

signal(SIGTERM, lambda signum, stack_frame: sys.exit(1))

def main():
    usage="""swga.py command
    
Available commands:
\t filter_primers: removes invalid primers from input
\t fg_locations:   finds and stores primer binding locations in foreground genome
\t make_graph:     creates the primer compatibility graph before finding sets
\t find_sets:      find initial sets of compatible primers
\t process_sets:   do additional filtering on compatible primer sets

"""
    parser = argparse.ArgumentParser(usage=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    parser.add_argument('-c', '--config_file', help='config file',
                        metavar="FILE")
                           
    args, remaining = parser.parse_known_args()
    fp_defaults = fl_defaults = mg_defaults = fs_defaults = ps_defaults = {}
    if args.config_file:
        config = ConfigParser.SafeConfigParser()
        config.read([args.config_file])
        fp_defaults = dict(config.items('filter_primers'))
        fl_defaults = dict(config.items('fg_locations'))
        mg_defaults = dict(config.items('make_graph'))
        fs_defaults = dict(config.items('find_sets'))
        ps_defaults = dict(config.items('process_sets'))

    subparsers = parser.add_subparsers(help='Commands', dest='command')

    filter_parser = subparsers.add_parser('filter_primers', help='''filter
    primers''', prog="filter_primers")
    fg_loc_parser = subparsers.add_parser('fg_locations', prog="fg_locations",
                                          help='''find primer binding locations
                                          in the foreground genome''') 
    mkgraph_parser = subparsers.add_parser('make_graph', prog='make_graph',
                                           help='''make initial primer
                                           compatibility graph''') 
    findsets_parser = subparsers.add_parser('find_sets', prog='find_sets',
                                            help='''find compatible sets of
                                            primers''') 
    processsets_parser = subparsers.add_parser('process_sets',
                                               prog='process_sets', help='''filter
                                               and analyze sets of primers''')



    # Filter primers command
    filter_parser.set_defaults(**fp_defaults)
    filter_parser.set_defaults(func=filter_primers)
    filter_parser.add_argument('--max_bg_binding', action='store',
                               type=int, help='''Max times a primer is allowed
                               to bind to the background genome.''')   
    filter_parser.add_argument('--num_primers', action='store', type=int,
                               help='''The max number of primers to use after
                               filtering.''')
    filter_parser.add_argument('primer_file', action='store',
                               help='''Space-delimited file where each row has
                               the primer sequence, fg binding number, bg
                               binding number, and fg/bg binding ratio, in that
                               order.''')


    # Foreground binding locations command
    fg_loc_parser.set_defaults(**fl_defaults)
    fg_loc_parser.set_defaults(func=fg_locations)
    fg_loc_parser.add_argument('--fg_genome', action='store',
                               help='''path to foreground genome
                               sequence, passed through the
                               utils/genome_flattener.sh program''')
    fg_loc_parser.add_argument('-i', '--input', action='store',
                               default=sys.stdin, type=argparse.FileType('r'),
                               help='''Space-delimited file where each row has
                               the primer sequence, fg binding number, bg
                               binding number, and fg/bg binding ratio, in that
                               order. If blank, reads from stdin.''', metavar='F')
    fg_loc_parser.add_argument('-o', '--output', action='store',
                               metavar='F', help='''Where to store the serialized output.''')
    fg_loc_parser.add_argument('--ncores', action='store', type=int, metavar='N',
                               default=multiprocessing.cpu_count(),
                               help='''Number of cores to use when
                               searching.''')
    fg_loc_parser.add_argument('-v', '--verbose', action='store_true',
                               help='''Display progress''')

    
    # Make primer graph command
    mkgraph_parser.set_defaults(**mg_defaults)
    mkgraph_parser.set_defaults(func=make_graph)
    mkgraph_parser.add_argument('--max_hetdimer_bind', type=int, help='''Max
    number of consecutive complimentary bases allowed between two primers.''')
    mkgraph_parser.add_argument('-i', '--input', action='store',
                                default=sys.stdin, type=argparse.FileType('r'),
                                help='''Space-delimited file where each row has
                                the primer sequence, fg binding number, bg
                                binding number, and fg/bg binding ratio, in that
                                order. If blank, reads from stdin.''')
    mkgraph_parser.add_argument('-o', '--output', action='store',
                                type=argparse.FileType('w', 0), 
                                default=sys.stdout, help='''Filename to store
                                the DIMACS-format output graph. If blank, writes
                                to stdout.''') 


    # Find sets command
    findsets_parser.set_defaults(**fs_defaults)
    findsets_parser.set_defaults(func=find_sets)
    findsets_parser.add_argument('-m', '--min_size', type=int,
                                 help='''Minimum size of primer sets.''')
    findsets_parser.add_argument('-M', '--max_size', type=int,
                                 help='''Maximum size of primer sets.''')
    findsets_parser.add_argument('-b', '--min_bg_bind_dist', type=int,
                                 help='''Minimum average distance between
                                 background binding sites for the primers in
                                 the set.''')
    findsets_parser.add_argument('-l', '--bg_genome_len', type=int,
                                 help='''Length of background genome.''')
    findsets_parser.add_argument('-s', '--set_finder',
                                 help='''Location of set_finder binary.''')

    findsets_parser.add_argument('-i', '--input', action='store',
                                 help='''Primer graph in DIMACS format where
                                 edges are between compatible primers. Default
                                 to stdin if unspecified.''',
                                 type=argparse.FileType('-r'),
                                 default=sys.stdin)
    findsets_parser.add_argument('-o', '--output', action='store',
                                 help='''Where to store results. Default to
                                 stdout if unspecified.''',
                                 type=argparse.FileType('w', 0),
                                 default=sys.stdout)


    # Filter sets command
    processsets_parser.set_defaults(**ps_defaults)
    processsets_parser.set_defaults(func=process_sets)
    processsets_parser.add_argument('-i', '--input', default=sys.stdin,
                                    type=argparse.FileType('r'),
                                    help='''Compatible sets of primers. One set
                                    per row, first number is the size of the
                                    set, following numbers are primer ids in
                                    that set, separated from spaces (output from
                                    find_sets command. Defaults to stdin if unspecified.''')  
    processsets_parser.add_argument('--fg_bind_locations',
                                    help='''Location of the sqlite database that
                                    contains foreground genome binding locations
                                    for each primer.''')  

    # parses the remaining subcommand options
    new_args = parser.parse_args(remaining)
    # calls the function corresponding to the subcommand with the specified options
    new_args.func(new_args)



def filter_primers(args):
    if not args.max_bg_binding:
        missing_default_value('max_bg_binding')
    elif not args.num_primers:
        missing_default_value('num_primers')
    filter_cmd = """sort -t ' ' -n -k 3 < {} | awk '{{if ($3 < {}) print
    $0}}' | sort -t ' ' -n -r -k 4 | head -n {}""".format(args.primer_file,
    args.max_bg_binding, args.num_primers) 
    subprocess.call(filter_cmd, shell=True, preexec_fn = lambda:
                        signal(SIGPIPE, SIG_DFL))

def fg_locations(args):
    primers = [line.strip('\n').split(' ')[0] for line in args.input]
    if args.verbose:
        print(("Populating genome binding locations for {} "
               "primers...").format(len(primers)))
    locations = ps.mp_find_primer_locations(primers, args.fg_genome,
                                            args.ncores, args.verbose)
    with gzip.GzipFile(args.output, 'w') as out:
        cPickle.dump(locations, out)
        if args.verbose:
            print("Locations serialized to {}".format(out.name))

        

def make_graph(args):
    pass
def find_sets(args):
    pass
def process_sets(args, parser):
    pass
    

        
def missing_default_value(missing_val):
    sys.stderr.write(("Error: {} not specified and no default found in config "
                     "file. Try -h for help.\n").format(missing_val))
    sys.exit(1)
                  
    


if __name__ == "__main__":
    main()
