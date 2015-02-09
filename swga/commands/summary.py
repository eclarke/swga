import os
import swga
import swga.database
from swga.commands import Command
from swga.database import Primer, Set
from swga.clint.textui import puts, indent, colored, max_width
from peewee import fn

def main(argv, cfg_file):
    cmd = Command('summary', cfg_file=cfg_file)
    cmd.parse_args(argv)
    summary(**cmd.args)


def summary(primer_db, fg_length, bg_length):
    
    db = swga.database.init_db(primer_db)
    db.connect()
    swga.database.create_tables(drop=False)

    avg_fg_bind, avg_bg_bind, nprimers = (Primer
                                          .select(fn.Avg(Primer.fg_freq),
                                                  fn.Avg(Primer.bg_freq),
                                                  fn.Count(Primer.seq))
                                          .scalar(as_tuple=True))
    fg_bind_ratio = avg_fg_bind / float(fg_length)
    bg_bind_ratio = avg_bg_bind / float(bg_length)
    nactive = Primer.select().where(Primer.active==True).count()
    min_tm, max_tm, avg_tm = (Primer
                              .select(fn.Min(Primer.tm),
                                      fn.Max(Primer.tm),
                                      fn.Avg(Primer.tm))
                              .where(Primer.active==True)
                              .scalar(as_tuple=True))

    nsets = Set.select(fn.Count(Set.sid)).scalar()
    if nsets > 0:
        bs = Set.select().order_by(Set.score.desc()).limit(1).get()
        bs_primers = ", ".join(swga.database.get_primers_for_set(bs.sid)).strip()
        best_set = bs.sid
        bs_size = bs.set_size
        bs_score = bs.score
        bs_stats = "- "+"\n - ".join("{}: {}".format(k,v)
                                     for k, v in bs.__dict__['_data'].items()
                             if k not in ["sid", "pids", "score"])

    summary_msg = """

    PRIMER SUMMARY
    ---------------
    There are {nprimers} primers in the database. 
    
    {nactive} are marked as active (i.e., they passed filter steps and will be used to find sets of compatible primers.) {ifzero_primers_msg}

    The average number of foreground genome binding sites is {avg_fg_bind:.0f}.
       (avg binding / genome_length = {fg_bind_ratio:05f})
    The average number of background genome binding sites is {avg_bg_bind:.0f}.
       (avg binding / genome_length = {bg_bind_ratio:05f})

    {melting_tmp_msg}


    SETS SUMMARY
    ---------------
    There are {nsets} sets in the database.
    {set_msg}---------------

    Report generated from {primer_db}
"""

    ifzero_primers_msg = colored.green("Run `swga filter` to identify primers to use." 
                                       if nactive == 0 else "")
    melting_tmp_msg =  ("""The melting temp of the primers ranges between {min_tm:.2f}C and {max_tm:.2f}C with an average of {avg_tm:.2f}C.
""" if nactive > 0 else "No melting temps have been calculated yet.")

    
    ifzero_sets_msg = colored.green("Run `swga sets` after identifying valid primers to "
                                    "begin collecting sets.\n")

    set_msg = ("""
    The best scoring set is #{best_set}, with {bs_size} primers and a score of {bs_score:03f}. Various statistics: 
    {bs_stats}
The primers in Set {best_set} are:
    {bs_primers}
    """ if nsets > 0 else ifzero_sets_msg)

    primer_db = os.path.abspath(primer_db)
    nprimers = colored.blue(nprimers, bold=True)
    nactive = colored.blue(nactive, bold=True)
    nsets = colored.blue(nsets, bold=True)
    set_msg = set_msg.format(**locals())
    melting_tmp_msg = melting_tmp_msg.format(**locals())
    summary_msg = summary_msg.format(**locals())

    with indent(2):
        puts(max_width(summary_msg, 80))
