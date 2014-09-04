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

    parser.add_argument('-i', '--input', 
                        default=sys.stdin, metavar='FILE',
                        type=argparse.FileType('r'), 
                        help="""Input file where each row contains a
                        primer, fg binding #, bg binding #, and fg/bg
                        binding ratio, separated by
                        whitespace. (default: stdin)""") 

    parser.add_argument('-g', '--genome', metavar='FILE',
                        required=True,
                        help="""Path to flattened genome
                        file. (default: %(default)s)""") 

    parser.add_argument('-l', '--locations_store', metavar='FILE',
                        help="""Where to store the primer binding
                        locations. (default: %(default)s)""") 

    parser.add_argument('-n', '--ncores', type=int, metavar="N",
                        default=multiprocessing.cpu_count(), 
                        help="""Number of cores to use when searching
                        (default: %(default)s).""") 

    parser.add_argument('-p', '--passthrough', action='store_true', 
                        help="""Echo input to stdout (for argument
                        chaining). (default: %(default)s)""") 

    args = parser.parse_args(argv)

    # Check to make sure genome file is valid
    if not os.path.isfile(args.genome):
        sys.stderr.write('Error: Genome specified by %s does not exist.\n' %
                         args.genome)
        exit(1)
    if not swga.check_if_flattened(args.genome):
        sys.stderr.write('Error: Genome does not appear to be flattened: use '
                         '`swga flatten` first.\n')
        exit(1)

    if not quiet and args.input.name == '<stdin>':
        swga.print_stdin_msg(parser.prog)

    if not quiet:
        swga.print_cfg_file(parser.prog, cfg_file)
        swga.print_args(parser.prog, args)

    primers = swga.read_primer_file(args.input, args.passthrough,
                                    not quiet)

    # Find locations using multiple processes
    locations = swga.mp_find_primer_locations(primers, args.genome,
                                              args.ncores, not quiet)

    # Save to gzipped pickled file (optimized for large numbers of sites)
    swga.save_locations(locations, args.locations_store, not quiet)


if __name__ == '__main__':
    main()
