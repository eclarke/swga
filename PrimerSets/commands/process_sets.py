import argparse
import ConfigParser
import os
import sys
import PrimerSets as ps
from functools import partial

def main(argv, cfg_file):
    defaults, _ = ps.parse_config(cfg_file, 'process_sets')
    parser = argparse.ArgumentParser(description="""Post-process and score sets
    outputted from find_sets.""")
    parser.set_defaults(**defaults)

    parser.add_argument('-i', '--input', default=sys.stdin,
    help='''Compatible sets of primers. One set per row, first number is the
    size of the set, following numbers are primer ids in that set, separated
    from spaces (output from find_sets command. (default: stdin)''')

    parser.add_argument('-o', '--output', default=sys.stdout,
                        type=argparse.FileType('w'), 
                        help='''Where to send output (default: stdout)''')

    parser.add_argument('--max_sets', type=int, help='''How many sets pass
    filter before we exit. (default: %(default)s)''')

    parser.add_argument('--max_fg_bind_dist', type=int, help='''Maximum distance
    between primers in a set on the foreground genome. (default: %(default)s)''')

    parser.add_argument('--fg_bind_locations', help='''Location of the output
    file that contains binding locations for each primer (from the fg_locations
    command). (default: %(default)s)''')

    score_funs = parser.add_mutually_exclusive_group()
    score_funs.add_argument('-s', '--score_expression', help="""Specify an
    expression to calculate the score added to the output. For variables
    available to use, see the README or docs. Incompatible with
    --plugin_score_fun argument. (default: %(default)s)""")

    score_funs.add_argument('-p', '--plugin_score_fun', help="""Specify a path
    to a function to use instead of the normal scoring function to create custom
    metrics or output. For help, see README or docs. Incompatible with
    --score_expression argument. (default: %(default)s)""")

    parser.add_argument('-v', '--verbose', action='store_true', help='''Display progress''')

    args = parser.parse_args(argv)
    process_sets(args)


def process_sets(args):
    '''
    Retrieves the primers and their binding locations from the output of
    find_sets and calculates the max binding distance between primers in the
    foreground genome.

    If the max distance is below a specified threshold, it passes the set
    and some additional attributes to a user-defined score function.

    After a specified number of sets pass the filter, it exits the process.
    '''
    primer_store = ps.load_locations(args.fg_bind_locations)
    # Find the user-defined scoring function
    score_fun = None
    if args.score_expression and args.plugin_score_fun:
        sys.stderr.write("Warning: User or config file specified both scoring "
        "expression and plugin score function. Using the plugin score function "
        "given by %s." % args.plugin_score_fun)
        score_fun = ps.get_user_fun(args.score_fun)
    elif args.score_expression:
        score_fun = partial(ps.default_score_set, expression=args.score_expression)
        # score_fun = lambda args: ps.default_score_set(args.score_expression, **args)

    passed = processed = 0
    for line in args.input:
        # Parse output from find_sets
        primer_ids, bg_ratio = ps.read_set_finder_line(line)
        primer_set = ps.get_primers_from_ids(primer_ids, primer_store)
        primer_locs = ps.get_primer_locations(primer_ids, primer_store)
        max_dist = max(ps.seq_diff(primer_locs))
        processed += 1
        if max_dist <= args.max_fg_bind_dist:
            passed += 1
            # Pass the set and attributes to the user-defined scoring function
            score_fun(primer_set=primer_set, primer_locs=primer_locs, max_dist=max_dist, bg_ratio=bg_ratio, output_handle=args.output)

        if args.verbose:
            sys.stderr.write("\rSets passing filter: \t{}/{}".format(passed, processed))
        if passed >= args.max_sets:
            sys.stderr.write("\nDone (scored %i sets). To quit, press Ctrl-C.")
            break
    if args.verbose:
        sys.stderr.write('\n')
    sys.exit()



if __name__ == '__main__':
    main()
