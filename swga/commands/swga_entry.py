import argparse

from swga import (
    error,
    DEFAULT_CFG_FNAME
)

from swga.commands import (
    init,
    summary,
    count,
    filter,
    find_sets,
    score,
    activate,
    export
)   

usage = """Usage: swga [-c --config CFG_FILE] <command> [options]

Utility commands:
  init:             initializes a directory with a pre-filled parameters file
  summary:          get a summary of the primers and sets found so far

Pipeline commands:
  count:            find kmer counts in foreground and background genomes
  filter:           filter kmers according to various criteria
  find_sets:        find and score compatible sets of primers in graph
  export:           export sets and primers

Other commands:
  activate:         activate a list of primers, calculating Tm and locations
  score:            score a set specified by a list of primers in a file

Options:
  --config FILE     path to config file (default %s)
"""


def main():
    command_opts = {
        'init': init.main,
        'summary': summary.main,
        'count': count.main,
        'filter': filter.main,
        'find_sets': find_sets.main,
        'score': score.main,
        'activate': activate.main,
        'export': export.main}

    cfg_file = DEFAULT_CFG_FNAME

    parser = argparse.ArgumentParser(
        usage=usage % cfg_file,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)

    parser.add_argument(
        'command',
        type=str,
        choices=command_opts.keys())

    parser.add_argument(
        '-c', '--config',
        metavar="CFG_FILE",
        help='path to config file (default: %(default)s)',
        default=cfg_file)

    args, remaining = parser.parse_known_args()

    try:
        command_opts[args.command](remaining, args.config)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)
