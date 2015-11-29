import functools
import subprocess
import os
import time
import signal

import click

import swga.graph as graph
import swga.locate as locate
import swga.score as score
import swga.sets as sets
from swga.primers import Primers
from swga import (warn, message)
from swga.commands._command import Command
from swga.workspace import Set

GRAPH_FP = "compatibility_graph.dimacs"
STATUS_LINE = '''\
\rSets: {: ^5,.6g} | Passed: {: ^5,.6g} | Smallest max bind dist:{: ^12,.4G}\
'''


class FindSets(Command):

    def run(self):
        if self.max_sets < 1:
            self.max_sets = float("inf")
        # We need to clear all the previously-used sets each time due to
        # uniqueness constraints
        all_sets = Set.select()
        if all_sets.count() > 0:
            if not self.force:
                click.confirm("Remove all previously-found sets?", abort=True)
            for s in all_sets:
                s.primers.clear()
                s.delete_instance()

        self.chr_ends = locate.chromosome_ends(self.fg_genome_fp)
        # Evaluate the scoring expression from a string and return it as a
        # callable function
        self.score_fun = functools.partial(
            score.default_score_set,
            expression=self.score_expression)

        graph.build_graph(self.max_dimer_bp, GRAPH_FP)

        message(
            "Finding sets. If nothing appears, try relaxing your parameters.")

        setfinder_lines = sets.find(
            min_bg_bind_dist=self.min_bg_bind_dist,
            min_size=self.min_size,
            max_size=self.max_size,
            bg_length=self.bg_length,
            graph_fp=GRAPH_FP,
            workers=self.workers)
        self.process_lines(setfinder_lines)

    def process_lines(self, setfinder_lines):
        passed = processed = 0
        smallest_max_dist = float('inf')

        try:
            for line in setfinder_lines:
                try:
                    primer_ids, bg_dist_mean = score.read_set_finder_line(line)
                except ValueError:
                    warn("Could not parse line:\n\t" + line)
                    continue

                primers = Primers.select_by_ids(primer_ids)
                processed += 1

                set_score, variables, max_dist = score.score_set(
                    primers=primers,
                    max_fg_bind_dist=self.max_fg_bind_dist,
                    bg_dist_mean=bg_dist_mean,
                    chr_ends=self.chr_ends,
                    score_fun=self.score_fun,
                    interactive=False
                )

                if max_dist < smallest_max_dist:
                    smallest_max_dist = max_dist

                message(
                    STATUS_LINE.format(processed, passed, smallest_max_dist),
                    newline=False)

                # Return early if the set doesn't pass
                if set_score is False:
                    continue
                else:
                    passed += 1

                Set.add(
                    _id=passed,
                    primers=primers,
                    score=set_score,
                    scoring_fn=self.score_expression,
                    **variables)

                if passed >= self.max_sets:
                    message("\nDone (scored %i sets)" % passed)
                    break
        finally:
            # Raises a GeneratorExit inside the find_sets command, prompting it
            # to quit the subprocess
            setfinder_lines.close()
