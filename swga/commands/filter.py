# -*- coding: utf-8 -*-
'''
Selects valid primers from a database according to various criteria.

This code could be a lot more succinct, but it's largely broken apart
to provide useful messages during operation (like how many primers satisfy
each criteria) and for performance reasons (i.e. we only calculate melt
temps on primers that pass the binding rate thresholds, and only calculate
binding locations on primers that pass all filters.)

Since all the results of the calculations are stored in the primer db,
the user can run this command multiple times to tune parameters and
subsequent runs will be much faster.
'''

import swga.primers
from swga.filters import Filter
from swga.commands import Command
from swga.database import Primer


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)

    # If we have an input file, use that. Otherwise pull from db
    if cmd.input:
        with open(cmd.input, 'rb') as infile:
            primers = swga.primers.read_primer_list(
                infile,
                cmd.fg_genome_fp,
                cmd.bg_genome_fp)
    else:
        cmd.skip_filtering = False
        primers = Primer.select()

    # Undo all active marks, if any
    deactivate_all_primers()

    if not cmd.skip_filtering:
        primers = (Filter(primers)
            .min_fg_rate(cmd.min_fg_bind)
            .max_bg_rate(cmd.max_bg_bind)
            .summarize()
            .tm_range(cmd.min_tm, cmd.max_tm)
            .limit_ratio(cmd.max_primers)
            .max_gini(cmd.max_gini, cmd.fg_genome_fp)
            .primers
        )

    n_active = activate_primers(primers)

    if n_active < cmd.max_primers:
        swga.warn(
            "Fewer than {} primers were selected ({} passed all the filters). "
            "You may want to try less restrictive filtering parameters."
            .format(cmd.max_primers, n_active))


def deactivate_all_primers():
    """Resets all active marks on primers."""
    Primer.update(active=False).execute()


def activate_primers(primers):
    """
    Marks as active all the kmers passed to it.
    """
    n_active = swga.primers.activate(primers)
    swga.message("Marked {} primers as active.".format(n_active))
    return n_active
