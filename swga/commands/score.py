import click
import functools
from peewee import fn
from swga import error, message
import swga.utils as utils
import swga.score as score
import swga.database
from swga.primers import Primers
from swga.commands._command import Command
from swga.database import Set


class Score(Command):

    def run(self):
        self.chr_ends = swga.locate.chromosome_ends(self.fg_genome_fp)
        # Evaluate the scoring expression from a string and return it as a
        # callable function
        self.score_fun = functools.partial(
            swga.score.default_score_set,
            expression=self.score_expression)

        primers = Primers(self.input)
        if len(primers) == 0:
            error(
                "No primers specified exist in database, aborting.",
                exception=False)

        bg_dist_mean = score.calculate_bg_dist_mean(primers, self.bg_length)

        set_score, variables, _ = score.score_set(
            primers=primers,
            max_fg_bind_dist=0,
            bg_dist_mean=bg_dist_mean,
            chr_ends=self.chr_ends,
            score_fun=self.score_fun,
            interactive=True
        )

        do_add_set, set_id = self.user_add_set(set_score, variables)

        if do_add_set:
            s = swga.database.add_set(
                _id=set_id,
                primers=primers,
                score=set_score,
                scoring_fn=self.score_expression,
                **variables)
            set_added = s is not None

            if set_added:
                swga.message("Set {} added successfully.".format(set_id))
            else:
                swga.message("That primer set already exists.")

    def user_add_set(self, set_score, variables):
        """Output set statistics and prompt the user to add the set."""
        set_dict = dict(
            {'score': set_score,
             'scoring_fn': self.score_expression}.items() +
            variables.items())
        message("Set statistics:\n - " + "\n - ".join(
            utils.fmtkv(k, v) for k, v in set_dict.items()))
        if self.force or (not self.force and click.confirm("Add set to database?",
                                                           default=True)):
            # User-provided sets have negative numbers, so we find the
            # smallest and decrement by 1
            min_set_id = Set.select(fn.Min(Set._id)).scalar()
            # This is None if there are no other sets yet
            if min_set_id is None:
                min_set_id = 0
            set_id = min_set_id - 1
            add_set = True
        else:
            add_set = False
        return add_set, set_id
