import functools
import subprocess
import os
import time
import signal

import click
import swga.database as database
import swga.graph as graph
import swga.locate as locate
import swga.score as score
import swga.utils as utils
from swga import (warn, message)
from swga.commands._command import Command
from swga.database import Set


GRAPH_FP = "compatibility_graph.dimacs"
STATUS_LINE = '''\
\rSets: {: ^5,.6g} | Passed: {: ^5,.6g} | Smallest max bind dist:{: ^12,.4G}\
'''


class FindSets(Command):

    def run(self):
        if self.max_sets < 1:
            self.max_sets = float("inf")
        # self.score_cmd = Score(argv)
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

        setfinder_lines = self.find_sets(GRAPH_FP, self.workers)
        self.process_lines(setfinder_lines)

    def find_sets(self, graph_fp, workers):
        if workers <= 1:
            return self._find_sets(graph_fp)
        else:
            return self._mp_find_sets(graph_fp, workers)

    def _find_sets(self, graph_fp, vertex_ordering="weighted-coloring"):
        assert vertex_ordering in ["weighted-coloring", "random"]
        find_set_cmd = [
            utils.set_finder, '-q', '-q',
            '--bg_freq', self.min_bg_bind_dist,
            '--bg_len', self.bg_length,
            '--min', self.min_size,
            '--max', self.max_size,
            '--unweighted',
            '--all',
            '--reorder', vertex_ordering,
            graph_fp
        ]
        find_set_cmd = " ".join([str(_) for _ in find_set_cmd])

        # We call the set_finder command as a line-buffered subprocess that
        # passes its output back to this process.
        # The function then yields each line as a generator; when close() is
        # called, it terminates the set_finder subprocess.
        process = subprocess.Popen(
            find_set_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            preexec_fn=os.setsid,
            bufsize=1)
        try:
            for line in iter(process.stdout.readline, b''):
                (yield line)
        finally:
            time.sleep(0.1)
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)

    def _mp_find_sets(self, graph_fp, workers):
        setfinder_procs = [self._find_sets(graph_fp, vertex_ordering="random")
                           for _ in range(workers)]
        try:
            while True:
                for i, setfinder in enumerate(setfinder_procs):
                    (yield setfinder.next())
        finally:
            for i, setfinder in enumerate(setfinder_procs):
                setfinder.close()

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

                primers = database.get_primers_for_ids(primer_ids)
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

                database.add_set(
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
