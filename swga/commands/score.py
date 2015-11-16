import click
import functools
from peewee import fn
from swga import error
import swga.locate as locate
import swga.score as score
import swga.database
from swga.primers import Primers
from swga.commands._command import Command
from swga.commands.summary import fmtkv
from swga.database import Set


class Score(Command):

    def __init__(self, argv):
        super(Score, self).__init__('score')
        self.parse_args(argv)
        self.chr_ends = swga.locate.chromosome_ends(self.fg_genome_fp)
        # Evaluate the scoring expression from a string and return it as a
        # callable function
        self.score_fun = functools.partial(
            swga.score.default_score_set,
            expression=self.score_expression)

    def run(self):
        primers = Primers(self.input)
        if len(primers) == 0:
            error(
                "No primers specified exist in database, aborting.",
                exception=False)

        self.score_set(
            set_id=0,
            primers=primers,
            # This is only used in filtering non-passing sets. If we're 
            # scoring manually, then this isn't relevant.
            max_fg_bind_dist=0,
            # bg_dist_mean is provided by the set_finder command. If we're
            # scoring a set manually, that's not supplied, so we calculate
            # it at runtime.
            bg_dist_mean=None,
            interactive=True, force=self.force)

    def score_set(self, set_id, primers, max_fg_bind_dist, bg_dist_mean,
                  interactive=False, force=True):

        binding_locations = locate.linearize_binding_sites(
            primers, self.chr_ends)
        max_dist = max(score.seq_diff(binding_locations))

        # If it's not a user-supplied set and it's not passing the filter,
        # abort immediately
        if not interactive and max_dist > max_fg_bind_dist:
            return False, max_dist

        if not bg_dist_mean:
            total_bg_freq = sum(p.bg_freq for p in primers)
            if total_bg_freq == 0:
                swga.warn(
                    "No primers appear in the background genome: "
                    "bg_dist_mean set as infinite")
                bg_dist_mean = float('Inf')
            else:
                bg_dist_mean = float(
                    self.bg_genome_len) / sum(p.bg_freq for p in primers)

        set_score, variables = self.score_fun(
            primer_set=primers,
            primer_locs=binding_locations,
            max_dist=max_dist,
            bg_dist_mean=bg_dist_mean)

        # Prompt the user if they want to add the set (if interactive)
        add_set = True
        if interactive:
            set_dict = dict(
                {'score': set_score,
                 'scoring_fn': self.score_expression}.items() +
                variables.items())
            swga.message("Set statistics:\n - " + "\n - ".join(
                fmtkv(k, v) for k, v in set_dict.items()))

            if force or (not force and click.confirm("Add set to database?",
                                                     default=True)):
                # User-provided sets have negative numbers, so we find the
                # smallest and decrement by 1
                min_set_id = Set.select(fn.Min(Set._id)).scalar()
                # This is None if there are no other sets yet
                if min_set_id is None:
                    min_set_id = 0
                self.set_id = min_set_id - 1
            else:
                add_set = False

        if add_set:
            s = swga.database.add_set(
                _id=set_id,
                primers=primers,
                score=set_score,
                scoring_fn=self.score_expression,
                **variables)
            set_added = s is not None

            if interactive and set_added:
                swga.message("Set {} added successfully.".format(set_id))
            elif interactive:
                swga.message("That primer set already exists.")

        return set_added, max_dist
