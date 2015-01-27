import sys
import swga
import swga.genome as genome
import swga.score as score
from swga.commands2 import Command
from functools import partial


def main(argv, cfg_file):
    cmd = Command('score', cfg_file=cfg_file)
    cmd.parse_args(argv)
    score_sets(**cmd.args)


def score_sets(input, 
               output,
               fg_bind_locations,
               score_expression,
               max_fg_bind_dist,
               max_sets,
               plugin_score_fun=None):
    '''
    Retrieves the primers and their binding locations from the output of
    find_sets and calculates the max binding distance between primers in the
    foreground genome.

    If the max distance is below a specified threshold, it passes the set
    and some additional attributes to a user-defined score function.

    After a specified number of sets pass the filter, it exits the process.
    '''
    primer_store = genome.load_locations(fg_bind_locations)
    # Find the user-defined scoring function
    score_fun = None
    if score_expression and plugin_score_fun:
        sys.stderr.write("Warning: User or config file specified both scoring "
                         "expression and plugin score function. Using "
                         "the plugin score function "
                         "given by %s." % plugin_score_fun)
        score_fun = swga.get_user_fun(score_fun)
    elif score_expression:
        score_fun = partial(score.default_score_set,
                            expression=score_expression)

    passed = processed = 0
    for line in input:
        # Parse output from find_sets
        primer_ids, bg_ratio = score.read_set_finder_line(line)
        primer_set = score.get_primers_from_ids(primer_ids, primer_store)
        primer_locs = score.get_primer_locations(primer_ids, primer_store)
        max_dist = max(score.seq_diff(primer_locs))
        processed += 1
        if max_dist <= max_fg_bind_dist:
            passed += 1
            # Pass the set and attributes to the user-defined scoring function
            score_fun(primer_set=primer_set,
                      primer_locs=primer_locs,
                      max_dist=max_dist,
                      bg_ratio=bg_ratio,
                      output_handle=output)

        sys.stderr.write(
            "\rSets passing filter: \t{}/{}".format(passed, processed))
        if passed >= max_sets:
            sys.stderr.write("\nDone (scored %i sets). To quit, press Ctrl-C.\n")
            break
    sys.exit()
