import os
import swga
from swga.commands2 import Command
from swga.primers import Primer, db, db_fname


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_kmers(**cmd.args)


def filter_kmers(kmer_dir,
                 fg_length,
                 bg_length,
                 min_avg_fg_bind,
                 max_avg_bg_bind,
                 max_bg_bind,
                 min_tm,
                 max_tm,
                 num_primers):

    primerdb_fp = os.path.join(kmer_dir, db_fname)
    
    if not os.path.isfile(primerdb_fp):
        swga.swga_error("Missing primers database: re-run `swga count`")

    fg_min_freq = float(min_avg_fg_bind) * fg_length
    bg_max_freq = float(max_avg_bg_bind) * bg_length

    # Only test for the higher of the two thresholds
    bg_max_freq = bg_max_freq if bg_max_freq > max_bg_bind else float(max_bg_bind)

    assert fg_min_freq > 0
    
    db.init(primerdb_fp)

    Primer.update(active = False).execute()
    
    valid_primers = (Primer
                     .update(active=True)
                     .where((Primer.fg_freq >= fg_min_freq) &
                            (Primer.bg_freq <= bg_max_freq) &
                            (Primer.tm <= max_tm) &
                            (Primer.tm >= min_tm))
                     .execute())
    
    swga.message("Marked {} qualifying primers".format(valid_primers))

    subquery = (Primer
                .select(Primer.seq)
                .where(Primer.active==True)
                .order_by(Primer.bg_freq)
                .limit(num_primers))

    query = (Primer
             .select()
             .where(Primer.seq << subquery)  # << means "IN"
             .order_by(Primer.ratio.desc()))
    
    for primer in query.execute():
        print primer.seq, primer.bg_freq, primer.ratio
