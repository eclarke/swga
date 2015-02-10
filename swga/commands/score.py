import sys
import json

import swga
import swga.primers
import swga.database
import swga.score as score
import swga.locate as locate
from swga.commands import Command
from functools import partial


def main(argv, cfg_file):
    cmd = Command('score', cfg_file=cfg_file)
    cmd.parse_args(argv)
    score_sets(**cmd.args)


def score_sets(primer_db, 
               fg_genome_fp,
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

    if max_sets < 1:
        max_sets = float("inf")

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

    swga.database.init_db(primer_db)

    chr_ends = locate.chromosome_ends(fg_genome_fp)
    passed = processed = 0
    for line in sys.stdin:
        try:
            primer_ids, bg_ratio = score.read_set_finder_line(line)
        except ValueError:
            swga.warn("Could not parse line:\n\t"+line)
            continue
        primers = swga.database.get_primers_for_ids(primer_ids)
        binding_locations = score.aggregate_primer_locations(primers) + chr_ends
        max_dist = max(score.seq_diff(binding_locations))
        processed += 1

        if max_dist <= max_fg_bind_dist:
            passed += 1
            set_score, variables = score_fun(primer_set=primers,
                                             primer_locs=binding_locations,
                                             max_dist=max_dist,
                                             bg_ratio=bg_ratio,
                                             output_handle=sys.stdout)
            swga.database.add_set(primers, score=set_score, 
                                 scoring_fn=score_expression, 
                                 pids=json.dumps(sorted(primer_ids)),
                                 **variables)
            
        sys.stderr.write(
            "\rSets passing filter: \t{}/{}".format(passed, processed))
        if passed >= max_sets:
            sys.stderr.write("\nDone (scored %i sets). To quit, press Ctrl-C.\n" % passed)
            sys.exit()
            break

