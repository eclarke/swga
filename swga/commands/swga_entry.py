import argparse
import swga.commands as commands
from swga import default_config_file, swga_error

usage="""Usage: swga [-c --config CFG_FILE] <command> [options]

Commands:
  init:             initializes a directory with a pre-filled parameters file
  autopilot:        runs a complete workflow based on parameters file
  summary:          get a summary of the primers and sets found so far

Other commands:
  count:            find kmer counts in foreground and background genomes
  filter:           filter kmers according to various criteria
  sets:             find compatible sets of primers in graph
  score:            score sets of primers
  export:           export information about primer sets in various formats

Options:
  --config FILE     path to config file (default %s)
  -q, --quiet       suppress messages
"""

def main():
    command_opts = {'init':commands.init.main,
                    'autopilot':commands.autopilot.main,
                    'summary':commands.summary.main,
                    'count':commands.count.main,
                    'filter':commands.filter.main,
                    'find_sets':commands.find_sets.main,
                    'score':commands.score.main,
                    'setopt':commands.setopt.main}

    cfg_file = default_config_file
    parser = argparse.ArgumentParser(usage=usage % cfg_file,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    parser.add_argument('command', type=str, choices=command_opts.keys())
    parser.add_argument('-c', '--config', metavar="CFG_FILE",
                        help='path to config file (default: %(default)s)',
                        default=cfg_file)
    args, remaining = parser.parse_known_args()
    try:
        command_opts[args.command](remaining, args.config)
    except KeyboardInterrupt:
        swga_error("\n-- Stopped by user --")


if __name__ == '__main__':
    main()
