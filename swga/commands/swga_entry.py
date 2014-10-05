import os
import sys
import argparse
import swga.commands as commands
from swga import default_config_file

usage="""Usage: swga [-c --config CFG_FILE] <command> [options]

Commands:
  count:            find kmer counts in foreground and background genomes
  filter:           filter kmers according to various criteria
  locate:           locate primers in a genome and store locations in file
  flatten:          flatten a genome file for easy searching
  mkgraph:          create heterodimer compatibility graph
  sets:             find compatible sets of primers in graph
  score:            score sets of primers
  locate-sets:      return locations of a set of primers in a genome
  autopilot:        run pipeline according to parameters given in config file

Options:
  --config FILE     path to config file (default %s)
  -q, --quiet       suppress messages
"""

def main():
    command_opts = {'count':commands.mer_count.main,
                    'flatten':commands.flatten.main,
                    'filter':commands.filter_primers.main,
                    'locate':commands.locate.main,
                    'mkgraph':commands.make_graph.main,
                    'sets':commands.find_sets.main,
                    'score':commands.process_sets.main,
                    'locate-sets':commands.locate_sets.main,
                    'autopilot':autopilot}
    cfg_file = default_config_file
    parser = argparse.ArgumentParser(usage=usage % cfg_file,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    parser.add_argument('command', type=str, choices=command_opts.keys())
    parser.add_argument('-c', '--config', metavar="CFG_FILE",
                        help='path to config file (default: %(default)s)',
                        default=cfg_file)
    parser.add_argument('-q', '--quiet', action='store_true')                        
    args, remaining = parser.parse_known_args()
    command_opts[args.command](remaining, args.config, args.quiet)


def autopilot(args, cfg_file):
    if not os.path.isfile(cfg_file):
        sys.stderr.write("Abort: no config file specified or the "
                         "specified file does not exist. Specify a "
                         "valid config file with --config")
        exit(1)


if __name__ == '__main__':
    main()
