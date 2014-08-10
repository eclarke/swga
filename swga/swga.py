#!/usr/bin/env python

import primer_sets as ps
import shlex
import sys
import argparse
import ConfigParser
import subprocess
import cPickle
import gzip
import os
import multiprocessing
import numpy as np
from multiprocessing.pool import ThreadPool
from signal import signal, SIGPIPE, SIG_DFL, SIGTERM



def main():
    usage="""swga.py [-c CONFIG_FILE] command

Available commands:
\t filter_primers: removes invalid primers from input
\t fg_locations:   finds and stores primer binding locations in foreground genome
\t make_graph:     creates the primer compatibility graph before finding sets
\t find_sets:      find initial sets of compatible primers
\t process_sets:   do additional filtering on compatible primer sets

By default, settings for all commands loaded from '{config}'. Override
config file location by setting the 'swga_params' environment
variable, or override certain settings by specifying them as arguments.

""".format(config=ps.default_config_file)

    
    parser = argparse.ArgumentParser(usage=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
                           
    args, remaining = parser.parse_known_args()
    fp_defaults = fl_defaults = mg_defaults = fs_defaults = ps_defaults = {}

    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', ps.default_config_file)
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
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
    fg_loc_parser.add_argument('--no_passthrough', action='store_true',
                               help='''Prevent input from passing through''')
    fg_loc_parser.add_argument('--fg_bind_locations', action='store',
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
                                 to stdin if unspecified.''', default='-')
    findsets_parser.add_argument('-o', '--output', action='store',
                                 help='''Where to store results. Default to
                                 stdout if unspecified.''', default='')


    # Filter sets command
    processsets_parser.set_defaults(**ps_defaults)
    processsets_parser.set_defaults(func=process_sets)
    processsets_parser.add_argument('-i', '--input', default=sys.stdin,
                                    help='''Compatible sets of primers. One set
                                    per row, first number is the size of the
                                    set, following numbers are primer ids in
                                    that set, separated from spaces (output from
                                    find_sets command. Defaults to
                                    stdin if unspecified.''')
    processsets_parser.add_argument('-o', '--output',
                                    default=sys.stdout,
                                    type=argparse.FileType('w'),
                                    help='''Where to send output''')
    processsets_parser.add_argument('--ncores', type=int,
                                    default=multiprocessing.cpu_count(),
                                    help='''Number of cores to use to do
                                    distance calculations''')
    processsets_parser.add_argument('--max_sets', type=int,
                                    help='''How many sets pass filter before we exit''')
    processsets_parser.add_argument('--max_fg_bind_dist', type=int,
                                    help='''Maximum distance between primers in
                                    a set on the foreground genome.''')
    processsets_parser.add_argument('--fg_bind_locations',
                                    help='''Location of the output file that
                                    contains foreground genome binding locations
                                    for each primer (from the fg_locations command).''')  
    processsets_parser.add_argument('-q', '--quiet', action='store_true', help='''Suppress progress output''')
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
    primers = []
    if args.no_passthrough:
        primers = ps.read_primers(args.input)
    else:
        pinput = [line for line in args.input]
        primers = ps.read_primers(pinput)
        sys.stdout.writelines(pinput)

    if args.verbose:
        sys.stderr.write(("Populating genome binding locations for {} "
                          "primers...\n").format(len(primers)))
    locations = ps.mp_find_primer_locations(primers, args.fg_genome,
                                            args.ncores, args.verbose)
    with gzip.GzipFile(args.fg_bind_locations, 'w') as out:
        cPickle.dump(locations, out)
        if args.verbose:
            sys.stderr.write("Locations stored in {}\n".format(out.name))


def make_graph(args):
    primers = ps.read_primers(args.input)
    arcs = ps.test_pairs(primers, args.max_hetdimer_bind)
    ps.write_graph(primers, arcs, args.output)

    
def find_sets(args):
    kwargs = vars(args)
    kwargs['output'] = '> '+args.output if args.output else ''
    find_set_cmd = ("{set_finder} -q -q -B {min_bg_bind_dist} -L {bg_genome_len}"
    " -m {min_size} -M {max_size} -a -u -r unweighted-coloring"
    " {input} {output}").format(**kwargs)
    subprocess.call(find_set_cmd, shell=True)
    

def process_sets(args):
    primer_locations = None
    passed = 0
    processed = 0
    with gzip.GzipFile(args.fg_bind_locations, 'r') as infile:
        primer_locations = cPickle.load(infile)
    for line in args.input:
        primer_set, primers, max_dist, stdev = ps.fg_bind_distances(line, primer_locations)
        primer_str = " ".join(primers)
        processed += 1
        if max_dist <= args.max_fg_bind_dist:
            passed += 1
            args.output.write("{} {} {}\n".format(stdev, max_dist,
                                                   primer_str))
        if not args.quiet:
            sys.stderr.write('\rSets passing filter: \t{}/{}'.format(passed, processed))
        if passed >= args.max_sets:
            sys.stderr.write('\nDone. If process fails to exit, press Ctrl-C. All data saved.')
            break
    if not args.quiet:
        sys.stderr.write('\n')
    sys.exit()

        
def missing_default_value(missing_val):
    sys.stderr.write(("Error: {} not specified and no default found in config "
                     "file. Try -h for help.\n").format(missing_val))
    sys.exit(1)
                  

if __name__ == "__main__":
    main()
