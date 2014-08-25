import argparse
import ConfigParser
import os
import sys
import stats
import gzip
import cPickle
import multiprocessing
import PrimerSets as ps


def main():
    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', ps.default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('process_sets'))

    parser = argparse.ArgumentParser(description="""Post-process and score sets
    outputted from find_sets.""")
    parser.set_defaults(**defaults)

    parser.add_argument('-i', '--input', default=sys.stdin,
    help='''Compatible sets of primers. One set per row, first number is the
    size of the set, following numbers are primer ids in that set, separated
    from spaces (output from find_sets command. (default: stdin)''')

    parser.add_argument('-o', '--output', default=sys.stdout,
    type=argparse.FileType('w'), help='''Where to send output (default: stdout)''')

    parser.add_argument('--max_sets', type=int, help='''How many sets pass
    filter before we exit''')

    parser.add_argument('--max_fg_bind_dist', type=int, help='''Maximum distance
    between primers in a set on the foreground genome.''')

    parser.add_argument('--fg_bind_locations', help='''Location of the output
    file that contains foreground genome binding locations for each primer
    (from the fg_locations command).''')

    parser.add_argument('-q', '--quiet', action='store_true', help='''Suppress
    progress output''')

    args = parser.parse_args()
    process_sets(args)


def process_sets(args):
    '''
    Takes the output from the set_finder binary and scores the sets,
    outputting only sets that pass a max foreground genome binding
    distance filter, specified in args.

    args: Namespace object from the argument parser
    '''
    primer_locations = None
    passed = 0
    processed = 0
    with gzip.GzipFile(args.fg_bind_locations, 'r') as infile:
        primer_locations = cPickle.load(infile)
    for line in args.input:
        primer_set, primers, max_dist, stdev = ps.fg_bind_distances(line,
        primer_locations, stats.stdev)
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

if __name__ == '__main__':
    main()
