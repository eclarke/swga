import click
import functools
from peewee import fn

import swga.primers
import swga.locate
import swga.score
import swga.database
from swga.primers import Primers
from swga.commands import Command
from swga.commands.summary import fmtkv
from swga.database import Set


def main(argv, cfg_file):
    cmd = Command('score')
    cmd.parse_args(argv)
    swga.database.init_db(cmd.primer_db)
    primers = Primers(cmd.input)
    if len(primers) == 0:
        swga.error("No primers specified exist in database, aborting.",
                   exception=False)
    chr_ends = swga.locate.chromosome_ends(cmd.fg_genome_fp)

    # Evaluate the user-defined scoring function
    score_fun = functools.partial(
        swga.score.default_score_set,
        expression=cmd.score_expression)

    score_set(
        set_id=0,
        primers=primers,
        chr_ends=chr_ends,
        score_fun=score_fun,
        score_expression=cmd.score_expression,
        bg_genome_len=cmd.bg_genome_len,
        bg_dist_mean=None,
        max_fg_bind_dist=0,
        interactive=True,
        force=cmd.force)


def score_set(
        set_id,
        bg_dist_mean,
        primers,
        chr_ends,
        score_fun,
        score_expression,
        max_fg_bind_dist,
        bg_genome_len=None,
        interactive=False,
        force=False):

    binding_locations = swga.locate.linearize_binding_sites(primers, chr_ends)
    max_dist = max(swga.score.seq_diff(binding_locations))

    # Abort now if it's not passing filter (and it's not a user-supplied set)
    if not interactive and max_dist > max_fg_bind_dist:
        return False, max_dist

    if not bg_dist_mean and not bg_genome_len:
        swga.error("Neither background length nor ratio were provided, "
                   "cannot calculate bg_dist_mean")
    elif not bg_dist_mean:
        total_bg_freq = sum(p.bg_freq for p in primers)
        if total_bg_freq == 0:
            swga.warn(
                "No primers appear in the background genome: bg_dist_mean set "
                "as infinite")
            bg_dist_mean = float('Inf')
        else:
            bg_dist_mean = float(
                bg_genome_len) / sum(p.bg_freq for p in primers)

    set_score, variables = score_fun(
        primer_set=primers,
        primer_locs=binding_locations,
        max_dist=max_dist,
        bg_dist_mean=bg_dist_mean)

    add_set = True
    # If it's user-supplied, they have the option of not adding it to db
    if interactive:
        set_dict = dict(
            {'score': set_score,
             'scoring_fn': score_expression}.items() + variables.items())
        swga.message("Set statistics:\n - " + "\n - ".join(
            fmtkv(k, v) for k, v in set_dict.items()))

        if force or (not force and click.confirm("Add set to database?", default=True)):
            # User-provided sets have negative numbers, so we find the smallest
            # and decrement by 1
            min_set_id = Set.select(fn.Min(Set._id)).scalar()
            # This is None if there are no other sets yet
            if min_set_id is None:
                min_set_id = 0
            set_id = min_set_id - 1
        else:
            add_set = False

    if add_set:
        s = swga.database.add_set(
            _id=set_id,
            primers=primers,
            score=set_score,
            scoring_fn=score_expression,
            **variables)
        set_added = s is not None

        if interactive and set_added:
            swga.message("Set {} added successfully.".format(set_id))
        elif interactive:
            swga.message("That primer set already exists.")
    return set_added, max_dist
