import argparse
import ConfigParser
import os
import re
import sys
from operator import itemgetter
from swga import default_config_file

def main():
    config = ConfigParser.SafeConfigParser()
    cfg_file = os.environ.get('swga_params', default_config_file)
    defaults = {}
    if os.path.isfile(cfg_file):
        config.read([cfg_file])
        defaults = dict(config.items('filter_primers'))

    parser = argparse.ArgumentParser(description="""Filter primers according to
    specified criteria.""")
    parser.set_defaults(**defaults)

    parser.add_argument('-M', '--max_bg_binding', action='store', type=int,
    help="""Max times a primer can bind to bg genome.""")

    parser.add_argument('-n', '--num_primers', action='store', type=int,
    help="""Max number of primers to use after filtering.""")

    parser.add_argument('-i', '--input', action='store',
    help="""Input file where each row contains a primer, fg binding #,
    bg binding #, and fg/bg binding ratio, separated by whitespace.""",
    default=sys.stdin, type=argparse.FileType('r'))

    parser.add_argument('-o', '--output', action='store',
    type=argparse.FileType('w',0), default=sys.stdout,
    help="""Filename to store the filtered primers (tab-delimited).""")

    parser.add_argument('-q', '--quiet', action='store_true',
    help="Suppress warnings")

    args = parser.parse_args()
    if not args.quiet and args.input.name == '<stdin>':
        sys.stderr.write("Receiving input from stdin...\n")
    primers = filter_primers(args)
    for line in primers:
        args.output.write("\t".join(line)+'\n')


def filter_primers(args):
    primers = []
    for i, line in enumerate(args.input):
        try:
            # use re.split so we can handle both tabs and spaces
            primer, fg_bind, bg_bind, ratio = re.split(r'[ \t]+', line.rstrip('\n'))
            primers.append([primer, int(fg_bind), int(bg_bind), float(ratio)])
        except ValueError as e:
            if not args.quiet:
                sys.stderr.write("Cannot parse line %i (reason: %s), skipping...\n" % (i, e.message))
            continue
    # Sort by bg binding numbers (want to select top n least bg binding)
    primers = sorted(primers, key=itemgetter(2))
    # Remove primers that bind too much to bg
    primers = [line for line in primers if line[2] <= args.max_bg_binding]
    # Sort by fg/bg ratio
    primers = sorted(primers, key=itemgetter(3))
    return primers[0:args.num_primers]


if __name__ == '__main__':
    main()
