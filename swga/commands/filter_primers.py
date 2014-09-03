import argparse
import sys
from operator import attrgetter
import PrimerSets as ps


def main(argv, cfg_file):
    defaults, _ = ps.parse_config(cfg_file, 'filter_primers')

    parser = argparse.ArgumentParser(description="""Filter primers according to
    specified criteria.""", prog='swga filter')
    parser.set_defaults(**defaults)

    parser.add_argument('-M', '--max_bg_binding', action='store', type=int,
    help="""Max times a primer can bind to bg genome. (default: %(default)s)""")

    parser.add_argument('-n', '--num_primers', action='store', type=int,
    help="""Max number of primers to use after filtering.
    (default: %(default)s)""")

    parser.add_argument('-i', '--input', action='store',
    help="""Input file where each row contains a primer, fg binding #,
    bg binding #, and fg/bg binding ratio, separated by whitespace.
    (default: stdin)""",
    default=sys.stdin, type=argparse.FileType('r'))

    parser.add_argument('-o', '--output', action='store',
    type=argparse.FileType('w',0), default=sys.stdout,
    help="""Filename to store the filtered primers (tab-delimited).
    (default: stdout)""")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Display messages")

    args = parser.parse_args(argv)
    if args.verbose and args.input.name == '<stdin>':
        sys.stderr.write("Receiving input from stdin...\n")
    primers = filter_primers(args)
    for primer in primers:
        args.output.write("{seq}\t{fg_freq}\t{bg_freq}\t{ratio}\n".format(**primer._asdict()))


def filter_primers(args):
    primers = ps.read_primer_file(args.input, False, args.verbose)
    # sort by bg binding count
    primers = sorted(primers, key=attrgetter("bg_freq"))
    # remove primers that bind too many times to bg
    primers = [p for p in primers if p.bg_freq <= args.max_bg_binding]
    # sort by fg/bg ratio
    primers = sorted(primers, key=attrgetter("ratio"))
    # return only the top <n>
    return primers[0:args.num_primers]


