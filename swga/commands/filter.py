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

import swga
import swga.primers
import swga.locate as locate
from swga.commands import Command
from swga.melting import Tm
from swga.database import Primer, update_in_chunks
from swga.clint.textui import progress

def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)
    
    # If we have an input file, use that. Otherwise pull from db
    if cmd.input:
        with open(cmd.input, 'rb') as infile:
            primers = swga.primers.read_primer_list(infile)
    else:
        primers = Primer.select()

    deactivate_all_primers()
    primers = filter_primers(primers, cmd.fg_min_avg_rate, cmd.bg_max_avg_rate,
                             cmd.fg_length, cmd.bg_length, cmd.min_tm,
                             cmd.max_tm, cmd.max_primers)
    update_locations(primers, cmd.fg_genome_fp)
    activate_primers(primers)


def deactivate_all_primers():
    """Resets all active marks on primers."""
    Primer.update(active=False).execute()

    
def activate_primers(primers):
    """
    Marks as active all the kmers passed to it.
    """
    n_active = Primer.update(active=True).where(Primer.seq << primers).execute()
    swga.message("Marked {} primers as active.".format(n_active))


def filter_primers(
        primers,
        fg_min_avg_rate,
        bg_max_avg_rate,
        fg_length,
        bg_length,
        min_tm,
        max_tm,
        max_primers):
    """
    Takes a list of sequences and retrieves them in the database, then returns
    those sequences that pass various criteria.
    """
    primers = Primer.select().where(Primer.seq << primers)
    fg_min_freq = float(fg_min_avg_rate) * fg_length
    bg_max_freq = float(bg_max_avg_rate) * bg_length

    # Find primers that pass the binding rate thresholds
    fgp = Primer.select().where((Primer.seq << primers) &
                                (Primer.fg_freq >= fg_min_freq))
    swga.message("{} primers bind foreground genome with avg rate >= {} sites/bp"
                 .format(fgp.count(), fg_min_avg_rate))
    
    bgp = Primer.select().where((Primer.seq << primers) &
                                (Primer.bg_freq <= bg_max_freq))
    swga.message("{} primers bind background genome with avg rate <= {} sites/bp"
                 .format(bgp.count(), bg_max_avg_rate))

    candidates = Primer.select().where((Primer.seq << primers) &
                                       (Primer.seq << fgp) &
                                       (Primer.seq << bgp))

    # Add melt temp for any primer that doesn't have it yet
    update_tms(candidates)
    
    valid_primers = Primer.select().where((Primer.seq << candidates) &
                                          (Primer.tm <= max_tm) &
                                          (Primer.tm >= min_tm))
    swga.message("{} of those primers have a melting temp within given range"
                 .format(candidates.count()))

    # Sort by background binding rate (smallest -> largest) and select top `n`,
    # then sort those by ratio (highest -> lowest)
    first_pass = (Primer.select()
                  .where(Primer.seq << valid_primers)
                  .order_by(Primer.bg_freq)
                  .limit(max_primers))

    second_pass = (Primer.select()
                   .where(Primer.seq << first_pass)
                   .order_by(Primer.ratio.desc()))

    return second_pass


def update_tms(primers):
    """
    Calculate the primer melting temps for any primers that don't already have
    them.
    Updates database with results.
    """
    primers = list(Primer.select().where((Primer.seq << primers) &
                                         (Primer.tm >> None)))
    if len(primers) > 0:
        for p in progress.bar(primers, label="Finding melting temps..."):
            p.tm = Tm(p.seq)
    update_in_chunks(primers, label="Updating database...")


def update_locations(primers, fg_genome_fp):
    """
    Find the primer binding sites in the foreground genome for primers that
    don't have the locations stored already.
    Updates database with results.
    """
    primers = list(Primer.select()
                   .where((Primer.seq << primers) &
                          (Primer.locations >> None))
                   .execute())
    # For more than 15 primers, we find the locations in parallel for performance
    if 0 < len(primers) < 15:
        for p in progress.bar(primers, label="Finding binding locations... "):
            p.locations = locate.primer_bind_sites(p.seq, fg_genome_fp)
    elif 0 < len(primers):
        primers = locate.primers_in_parallel(primers,
                                             fg_genome_fp)
    update_in_chunks(primers, label = "Updating database... ")
    
