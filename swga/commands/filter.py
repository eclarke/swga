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

from swga.primers import Primers
from swga.commands import Command
import swga.database


def main(argv, cfg_file):
    cmd = Command('filter')
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)

    # If we have an input file, use that. Otherwise pull from db
    if cmd.input:
        with open(cmd.input, 'rb') as infile:
            primers = Primers(infile)
    else:
        cmd.skip_filtering = False
        primers = Primers()

    assert isinstance(primers, Primers)

    # Undo all active marks, if any
    swga.database.Primer.update(active=False).execute()

    if not cmd.skip_filtering:
        (
            primers
            .filter_min_fg_rate(cmd.min_fg_bind)
            .filter_max_bg_rate(cmd.max_bg_bind)
            .summarize()
            .filter_tm_range(cmd.min_tm, cmd.max_tm)
            .limit_to(cmd.max_primers)
            .filter_max_gini(cmd.max_gini, cmd.fg_genome_fp)
        )

    primers.activate(cmd.max_primers)
