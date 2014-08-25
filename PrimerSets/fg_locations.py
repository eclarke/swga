import argparse
import ConfigParser
import os
import sys
import multiprocessing
import cPickle
import gzip
import PrimerSets

def main():
    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', PrimerSets.default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('fg_locations'))

    parser = argparse.ArgumentParser(description="""Find binding locations of
    primers in the foreground genome.""")
    parser.set_defaults(**defaults)

    parser.add_argument('-i', '--input', default=sys.stdin, metavar='FILE',
    type=argparse.FileType('r'), help="""Input file where each row contains a
    primer, fg binding #, bg binding #, and fg/bg binding ratio, separated by
    whitespace. (default: stdin)""")

    parser.add_argument('-f', '--fg_genome', metavar='FILE',
    help="""Path to flattened foreground genome file. (default: %(default)s)""")

    parser.add_argument('-l', '--fg_bind_locations', metavar='FILE',
    help="Where to store the binding locations. (default: %(default)s)")

    parser.add_argument('-n', '--ncores', type=int, metavar="N",
    default=multiprocessing.cpu_count(), help="""Number of cores to use when
    searching (default: %(default)s).""")

    parser.add_argument('-p', '--passthrough', help="""Echo input to stdout
    (for argument chaining). (default: %(default)s)""", action='store_true')

    parser.add_argument('-v', '--verbose', action='store_true', help="""Display
    progress.""")

    args = parser.parse_args()
    if args.verbose and args.input.name == '<stdin>':
        sys.stderr.write("Receiving input from stdin...\n")

    primers = PrimerSets.read_primer_file(args.input, args.passthrough,
    args.verbose)

    locations = PrimerSets.mp_find_primer_locations(primers, args.fg_genome,
    args.ncores, args.verbose)

    with gzip.GzipFile(args.fg_bind_locations, 'w') as out:
        cPickle.dump(locations, out)
        if args.verbose:
            sys.stderr.write("Locations stored in {}\n".format(out.name))



if __name__ == '__main__':
    main()
