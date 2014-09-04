import argparse
import ConfigParser
import os
import sys
import multiprocessing
import cPickle
import gzip
import swga


def main(argv, cfg_file, quiet):
    '''
    Find binding locations of primers in a genome and create a
    primer store at the specified location.
    '''
    parser = swga.basic_cmd_parser(description=main.__doc__, 
                                   cmd_name="locate",
                                   cfg_file=cfg_file)

    parser.add_argument('-i', '--input', default=sys.stdin, metavar='FILE',
    type=argparse.FileType('r'), help="""Input file where each row contains a
    primer, fg binding #, bg binding #, and fg/bg binding ratio, separated by
    whitespace. (default: stdin)""")

    parser.add_argument('-g', '--genome', metavar='FILE',
    help="""Path to flattened genome file. (default: %(default)s)""")

    parser.add_argument('-l', '--locations_store', metavar='FILE',
    help="Where to store the primer binding locations. (default: %(default)s)")

    parser.add_argument('-n', '--ncores', type=int, metavar="N",
    default=multiprocessing.cpu_count(), help="""Number of cores to use when
    searching (default: %(default)s).""")

    parser.add_argument('-p', '--passthrough', help="""Echo input to stdout
    (for argument chaining). (default: %(default)s)""", action='store_true')

    args = parser.parse_args(argv)
    if not quiet and args.input.name == '<stdin>':
        swga.print_stdin_msg(parser.prog)

    # Check to make sure foreground genome is valid
    if not os.path.isfile(args.fg_genome):
        raise ValueError('Genome specified by %s does not exist.' %
                         args.fg_genome)
    if not swga.check_if_flattened(args.fg_genome):
        raise ValueError('Genome does not appear to be flattened: use '
                         '`swga flatten` first')

    if not quiet and args.input.name == '<stdin>':
        sys.stderr.write("Receiving input from stdin...\n")

    primers = swga.read_primer_file(args.input, args.passthrough,
    not args.quiet)

    # Find locations using multiple processes
    locations = swga.mp_find_primer_locations(primers, args.fg_genome,
    args.ncores, not args.quiet)

    # Save to gzipped pickled file (optimized for large numbers of sites)
    swga.save_locations(locations, args.fg_bind_locations, not args.quiet)



if __name__ == '__main__':
    main()
