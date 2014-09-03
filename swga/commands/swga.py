import os
import sys
import argparse
import PrimerSets.commands as commands
import PrimerSets as ps

usage="""Usage: swga [-c --config CFG_FILE] <command> [options]

Commands:
  count:            find kmer counts in foreground and background genomes
  filter:           filter kmers according to various criteria
  flatten:          flatten a genome file for easy searching
  mkgraph:          create heterodimer compatibility graph
  sets:             find compatible sets of primers in graph
  score:            score sets of primers
  autopilot:        run pipeline according to parameters given in config file

Options:
  --config FILE     path to config file (default %s)
  -v, --verbose     display messages
"""

def main():
    command_opts = {'count':commands.count_mers.main,
                    'filter':commands.filter_primers.main,
                    'flatten':commands.flatten.main,
                    'mkgraph':commands.make_graph.main,
                    'sets':commands.find_sets.main,
                    'score':commands.process_sets.main,
                    'autopilot':autopilot}
    swgahome = ps.get_swgahome()
    if os.path.isfile(ps.default_config_file):
        cfg_file = ps.default_config_file
    else:
        cfg_file = os.path.join(swgahome, 'parameters.cfg')
    parser = argparse.ArgumentParser(usage=usage % cfg_file,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    parser.add_argument('command', type=str, choices=command_opts.keys())
    parser.add_argument('-c', '--config', metavar="CFG_FILE",
                        help='path to config file (default: %(default)s)',
                        default=cfg_file)
    parser.add_argument('-v', '--verbose', action='store_true')                        
    args, remaining = parser.parse_known_args()
    if args.verbose:
        sys.stderr.write("swga %s: using config file at %s\n" %
                         (args.command, os.path.abspath(args.config)))
    command_opts[args.command](remaining, args.config)


def autopilot(args, cfg_file):
    if not os.path.isfile(cfg_file):
        sys.stderr.write("Abort: no config file specified or the "
                         "specified file does not exist. Specify a "
                         "valid config file with --config")
        exit(1)


if __name__ == '__main__':
    main()
