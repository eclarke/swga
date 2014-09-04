import argparse
import sys
from operator import attrgetter
import swga


def main(argv, cfg_file, quiet):
    '''Filter primers according to specified criteria.'''
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='filter',
                                   cfg_file=cfg_file)

    parser.add_argument('-M', '--max_bg_binding', action='store', type=int,
                        help="""Max times a primer can bind to bg
                        genome. (default: %(default)s)""") 

    parser.add_argument('-n', '--num_primers', action='store', type=int,
                        help="""Max number of primers to use after
                        filtering. (default: %(default)s)""") 

    parser.add_argument('-i', '--input', action='store',
                        default=sys.stdin, type=argparse.FileType('r'),
                        help="""Input file where each row contains a
                        primer, fg binding #, bg binding #, and fg/bg
                        binding ratio, separated by
                        whitespace. (default: stdin)""")

    parser.add_argument('-o', '--output', action='store',
                        type=argparse.FileType('w',0), default=sys.stdout,
                        help="""Filename to store the filtered primers
                        (tab-delimited). (default: stdout)""") 

    args = parser.parse_args(argv)
    if not quiet and args.input.name == '<stdin>':
        swga.print_stdin_msg(parser.prog)
    
    primers = filter_primers(args, quiet)
    for primer in primers:
        args.output.write("{seq}\t{fg_freq}\t{bg_freq}\t{ratio}\n".format(**primer._asdict()))


def filter_primers(args, quiet):
    primers = swga.read_primer_file(args.input, False, not quiet)
    # sort by bg binding count
    primers = sorted(primers, key=attrgetter("bg_freq"))
    # remove primers that bind too many times to bg
    primers = [p for p in primers if p.bg_freq <= args.max_bg_binding]
    # keep only the top <num_primers>
    primers = primers[0:args.num_primers]
    # sort by fg/bg ratio
    primers = sorted(primers, key=attrgetter("ratio"))
    return primers
    


