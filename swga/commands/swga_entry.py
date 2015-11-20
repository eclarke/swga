import argparse
import os
from swga import (
    error,
    DEFAULT_DB_FNAME,
    DEFAULT_CFG_FNAME
)
import swga.database as database
from swga.commands import (
    init,
    Summary,
    Count,
    Filter,
    FindSets,
    Score,
    Activate,
    Export
)

usage = """Usage: swga <command> [options]

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
"""


def main():
    command_opts = {
        'init': init.main,
        'summary': Summary,
        'count': Count,
        'filter': Filter,
        'find_sets': FindSets,
        'score': Score,
        'activate': Activate,
        'export': Export
    }

    parser = argparse.ArgumentParser(
        usage=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)

    parser.add_argument(
        'command',
        type=str,
        choices=command_opts.keys())

    args, remaining = parser.parse_known_args()

    try:
        if args.command == 'init':
            init.main(remaining)
        else:
            db_name = os.path.abspath(DEFAULT_DB_FNAME)
            cfg_file = os.path.abspath(DEFAULT_CFG_FNAME)
            with database.connection(db_name):
                database.check_version()
                metadata = database.get_metadata(db_name)
                cmd = command_opts[args.command](args.command, cfg_file, metadata)
                cmd.parse_args(remaining)
                cmd.run()
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)
