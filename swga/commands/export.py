import argparse
import swga
from pyfasta import Fasta
from swga import genome
import sys
import os
from swga import stats
from collections import Counter

def main(argv, cfg_file, quiet):
    """
    Export the locations of a primer set in a target genome in various formats.
    """
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name="export",
                                   cfg_file=cfg_file)
    export_opts = {'bed':export_bed,
                   'bedgraph':export_bedgraph}
    parser.add_argument('-f', '--format', required=True,
                        help="One of {}".format(export_opts.keys()))
    parser.add_argument('-i', '--input', required=True,
                        default=sys.stdin, metavar="FILE",
                        type=argparse.FileType('r'),
                        help="""File with primer sequences, one per line 
                        (default: stdin)""")
    parser.add_argument('-t', '--target', required=True,
                        help="""Path to target genome in FASTA format (not 
                        flattened)""")
    parser.add_argument('-o', '--output_folder',
                        required=True,
                        help="Path to folder to store output files")
    parser.add_argument('-b', '--opts_str',
                        metavar="OPTS",
                        help="""(for bedgraph) Additional parameters for 
                        BedGraph file format, in quotes, that will be inserted 
                        after `track type=bedGraph` in the output file.""")
    parser.add_argument('-s', '--sliding_window_size',
                        metavar="SIZE",
                        type=int,
                        default=10000,
                        help="""(for bedgraph) The size of the sliding window 
                        used to calculate average primer binding coverage 
                        (default: %(default)s)""")
    parser.add_argument('-w', '--sliding_window_move',
                        metavar="BASES",
                        type=int,
                        help="""(for bedgraph) The number of bases by which the 
                        sliding window moves each time. (default: window_size/5)""")
    
    args = parser.parse_args(argv)
    primers = [_.strip('\n') for _ in args.input.readlines()]
    target = Fasta(args.target)
    export_opts[args.format](args, primers, target)

def export_bed(args, primers, target):
    target_name = '.'.join(args.target.split('.')[:-1])
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    with open(os.path.join(args.output_folder, 'whole_set.bed'), 'w') as whole_set_file:
        whole_set_file.write("track name={}-whole-set\n".format(args.output_folder))
        for primer in primers:
            fname = os.path.join(args.output_folder, 
                                 "{}-{}.bed".format(target_name, primer))
            with open(fname, 'w') as bedfile:
                bedfile.write("track name={}\n".format(primer))
                for record in target:
                    seq = target[record][::]
                    record_name = record.split('|')[0].strip()
                    for location in genome.find_locations(primer, seq):
                        record_string = "{} {} {}\n".format(record_name, 
                                                            location,
                                                            location+len(primer))
                        bedfile.write(record_string)
                        whole_set_file.write(record_string)


def export_bedgraph(args, primers, target):
    target_name = '.'.join(args.target.split('.')[:-1])
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    with open(os.path.join(args.output_folder, '{}.bedgraph'.format(target_name)), 'w') as outfile:
        typestr = "track type=bedGraph {}\n".format(args.opts_str)
        outfile.write(typestr)
        s = args.sliding_window_size
        move = args.sliding_window_move if args.sliding_window_move else int(s/5)
        for record in target:
            seq = target[record][::]
            record_name = record.split('|')[0].strip()
            # Counter objects tally the unique items assigned to them
            tmp_locs = Counter()
            for primer in primers:
                for l in genome.find_locations(primer, seq):
                    tmp_locs += Counter(xrange(l, l+len(primer)))                                        
            # Gives the number of times a primer binds to that site
            hits = [tmp_locs[i] for i in xrange(0, len(seq))]
            for start in xrange(0, len(seq)-s, move):
                end = start+s
                avg_binding = stats.mean(hits[start:end])
                linestr = "{} {} {} {}\n".format(record_name, start, end, avg_binding)
                outfile.write(linestr)
