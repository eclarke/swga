import os
import sys
import argparse
import PrimerSets.scripts as scripts
import PrimerSets as ps

usage="""Usage: swga [-c --config CFG_FILE] <command> [options]

Commands:
  count:            find kmer counts in foreground and background genomes
  filter:           filter kmers according to various criteria
  flatten:          flatten a genome file for easy searching
  sets:             find compatible sets of primers
  score:            score sets of primers
  autopilot:        run pipeline according to parameters given in config file

Options:
  --config FILE     path to config file (default %s)
"""

def main():
    commands = {
        'count':scripts.count_mers.main,
        'filter':scripts.filter_primers.main,
#        'flatten':scripts.flatten.main,
        'sets':scripts.find_sets.main,
        'score':scripts.process_sets.main,
        'autopilot':autopilot}
    cfg_file = os.environ.get('swga_params', ps.default_config_file)
    parser = argparse.ArgumentParser(
        usage=usage % cfg_file,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)
    parser.add_argument('command', type=str, choices=commands.keys())
    parser.add_argument('-c', '--config', metavar="CFG_FILE",
        help='path to config file (default: %(default)s)',
        default=cfg_file)
    args, remaining = parser.parse_known_args()
    commands[args.command](remaining, args.config)


def autopilot(args, cfg_file):
    if not os.path.isfile(cfg_file):
        sys.stderr.write("Abort: no config file specified or does not exist. "
        "Specify a config file path with --config or set the $swga_params "
        "environment variable.\n")
        sys.exit(1)





if __name__ == '__main__':
    main()
