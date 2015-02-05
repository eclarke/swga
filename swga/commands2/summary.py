import os
import swga
from swga.commands2 import Command
from swga.primers import Primer, Set
from swga.clint.textui import puts, indent, colored, max_width
import swga.stats as stats


def main(argv, cfg_file):
    cmd = Command('summary', cfg_file=cfg_file)
    cmd.parse_args(argv, quiet=True)
    summary(**cmd.args)


def summary(primer_db, fg_length, bg_length):
    
    db = swga.primers.init_db(primer_db)
    db.connect()
    swga.primers.create_tables(drop=False)

    primers = Primer.select()
    avg_fg_bind = stats.mean([p.fg_freq for p in primers])
    fg_bind_ratio = avg_fg_bind / float(fg_length)
    avg_bg_bind = stats.mean([p.bg_freq for p in primers])
    bg_bind_ratio = avg_bg_bind / float(bg_length)
    nactive = Primer.select().where(Primer.active==True).count()
    tms = [p.tm for p in primers]
    min_tm = min(tms)
    max_tm = max(tms)
    avg_tm = stats.mean(tms)


    nsets = Set.select().count()
    if nsets > 0:
        bs = Set.select().order_by(Set.score.desc()).limit(1).get()
        bs_primers = "\n    ".join(swga.primers.get_primers_for_set(bs.sid))
        best_set = bs.sid
        bs_size = bs.set_size
        bs_score = bs.score


    summary_msg = """

    PRIMER SUMMARY
    ---------------
    There are {nprimers} primers in the database. 
    
    {nactive} are marked as active (i.e., they passed filter steps and will be used to find sets of compatible primers.) {ifzero_primers_msg}

    The average number of foreground genome binding sites is {avg_fg_bind:.0f}.
       (avg binding / genome_length = {fg_bind_ratio:05f})
    The average number of background genome binding sites is {avg_bg_bind:.0f}.
       (avg binding / genome_length = {bg_bind_ratio:05f})

    The melting temp of the primers ranges between {min_tm:.2f}C and
    {max_tm:.2f}C with an average of {avg_tm:.2f}C.


    SETS SUMMARY
    ---------------
    There are {nsets} sets in the database.
    {set_msg}
    ---------------

    Report generated from {primer_db}
"""

    ifzero_primers_msg = colored.green("Run `swga filter` to identify primers to use." 
                                       if nactive == 0 else "")
    ifzero_sets_msg = colored.green("Run `swga sets` after identifying valid primers to "
                                    "begin collecting sets.")

    set_msg = ("""
    The best scoring set is #{best_set}, with {bs_size} primers and a score of {bs_score:.2f}. The primers in Set {best_set} are:
    {bs_primers}
    """ if nsets > 0 else ifzero_sets_msg)

    primer_db = os.path.abspath(primer_db)
    nprimers = colored.blue(len(tms), bold=True)
    nactive = colored.blue(nactive, bold=True)
    nsets = colored.blue(nsets, bold=True)
    set_msg = set_msg.format(**locals())
    summary_msg = summary_msg.format(**locals())

    with indent(2):
        puts(max_width(summary_msg, 80))
