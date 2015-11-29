import argparse
import os
from swga import (
    error,
    DEFAULT_DB_FNAME,
    DEFAULT_CFG_FNAME,
    __version__
)
import workspace
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


def setup_and_run(cmd_class, name, remaining_args):
    """Setup and run a command in the current workspace."""
    db_name = os.path.abspath(DEFAULT_DB_FNAME)
    cfg_file = os.path.abspath(DEFAULT_CFG_FNAME)
    with workspace.connection(db_name) as ws:
        assert not ws.is_closed()
        ws.check_version(__version__)
        metadata = ws.metadata
        cmd = cmd_class(name, cfg_file, metadata)
        cmd.parse_args(remaining_args)
        cmd.run()


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
            command_opts[args.command](remaining)
        else:
            cmd_class = command_opts[args.command]
            setup_and_run(cmd_class, args.command, remaining)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)
