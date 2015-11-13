import functools

import click
import swga.database
import swga.graph as graph
import swga.locate as locate
import swga.score as score
from swga.primers import Primers
from swga.clint.textui import progress
from swga.commands import Command
from swga.commands.score import score_set
from swga.database import Set, init_db
import swga.setfinder as setfinder

graph_fname = "compatibility_graph.dimacs"


def main(argv, cfg_file):
    cmd = Command('find_sets')
    score_cmd = Command('score')
    cmd.parse_args(argv)
    score_cmd.parse_args(argv)

    init_db(cmd.primer_db)

    # We need to clear all the previously-used sets each time due to uniqueness
    # constraints
    allsets = Set.select()
    if allsets.count() > 0:
        if not cmd.force:
            click.confirm("Remove all previously-found sets?", abort=True)
        for s in progress.bar(allsets, expected_size=allsets.count()):
            s.primers.clear()
            s.delete_instance()

    make_graph(cmd.max_dimer_bp, graph_fname)

    swga.message(
        "Now finding sets. If nothing appears, try relaxing your parameters."
    )

    if cmd.workers <= 1:
        setlines = setfinder.find_sets(
            cmd.min_bg_bind_dist,
            cmd.min_size,
            cmd.max_size,
            cmd.bg_genome_len,
            graph_fp=graph_fname)
    else:
        setlines = setfinder.mp_find_sets(
            nprocesses=cmd.workers,
            graph_fp=graph_fname,
            min_bg_bind_dist=cmd.min_bg_bind_dist,
            min_size=cmd.min_size,
            max_size=cmd.max_size,
            bg_genome_len=cmd.bg_genome_len)

    score_sets(
        setlines,
        cmd.fg_genome_fp,
        score_cmd.score_expression,
        cmd.max_fg_bind_dist,
        cmd.max_sets)


def make_graph(max_hetdimer_bind, outfile):
    '''Selects all active primers and outputs a primer compatibility graph.'''

    # Reset all the primer IDs (as ids are only used for set_finder)
    primers = Primers.select_active().assign_ids()
    print [(p._id, p.ratio) for p in primers]
    swga.message("Composing primer compatibility graph...")
    edges = graph.test_pairs(primers, max_hetdimer_bind)
    print edges

    if len(edges) == 0:
        swga.error(
            "No compatible primers. Try relaxing your parameters.",
            exception=False)

    with open(outfile, 'wb') as out:
        graph.write_graph(primers, edges, out)


def score_sets(
        setlines,
        fg_genome_fp,
        score_expression,
        max_fg_bind_dist,
        max_sets):
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

    # Evaluate the user-defined scoring function
    score_fun = functools.partial(
        score.default_score_set,
        expression=score_expression)

    chr_ends = locate.chromosome_ends(fg_genome_fp)

    passed = processed = 0
    min_max_dist = float('inf')

    status = "\rSets: {: ^5,.6g} | Passed: {: ^5,.6g} | Smallest max binding distance:{: ^12,.4G}"

    try:
        for line in setlines:
            try:
                primer_ids, bg_dist_mean = score.read_set_finder_line(line)
            except ValueError:
                swga.warn("Could not parse line:\n\t" + line)
                continue

            primers = swga.database.get_primers_for_ids(primer_ids)
            processed += 1

            set_passed, max_dist = score_set(
                set_id=passed,
                bg_dist_mean=bg_dist_mean,
                primers=primers,
                chr_ends=chr_ends,
                score_fun=score_fun,
                score_expression=score_expression,
                max_fg_bind_dist=max_fg_bind_dist,
                interactive=False)

            if set_passed:
                passed += 1

            min_max_dist = max_dist if max_dist < min_max_dist else min_max_dist

            swga.message(
                status.format(processed, passed, min_max_dist), newline=False)

            if passed >= max_sets:
                swga.message("\nDone (scored %i sets)" % passed)
                break
    finally:
        # Raises a GeneratorExit inside the find_sets command, prompting it to quit
        # the subprocess
        setlines.close()
