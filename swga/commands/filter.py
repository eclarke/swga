
import swga
import swga.locate as locate
from swga.core import chunk_iterator
from swga.commands import Command
from swga.melting import Tm
from swga.primers import Primer, upsert_chunk
from swga.clint.textui import progress

def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_kmers(**cmd.args)


def filter_kmers(primer_db,
                 fg_genome_fp,
                 fg_length,
                 bg_length,
                 min_avg_fg_rate,
                 max_avg_bg_rate,
                 max_bg_bind,
                 min_tm,
                 max_tm,
                 num_primers):
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
    
    swga.primers.init_db(primer_db)

    fg_min_freq = float(min_avg_fg_rate) * fg_length
    bg_max_freq = float(max_avg_bg_rate) * bg_length

    # An 'active' flag means the primer passed all the criteria. We reset this
    # each run in case parameters change.
    Primer.update(active = False).execute()

    # Find primers that pass binding rate thresholds (the "candidates")
    fgp = Primer.select(Primer.seq).where(Primer.fg_freq >= fg_min_freq)
    swga.message("{} primers bind foreground genome with avg rate >= {} sites/bp"
                 .format(fgp.count(), min_avg_fg_rate)) 

    bgp = Primer.select(Primer.seq).where(Primer.bg_freq <= bg_max_freq)
    swga.message("{} primers bind background genome with avg rate <= {} sites/bp"
                 .format(bgp.count(), max_avg_bg_rate))

    missing_tms = list(Primer
                      .select()
                      .where((Primer.seq << fgp) & (Primer.seq << bgp) &
                             (Primer.tm >> None)))

    # Calculate all the melting temps in one go, then update the db in chunks
    # for performance reasons
    if len(missing_tms) > 0:
        for p in progress.bar(missing_tms, label="Finding melting temps..."):
            p.tm = Tm(p.seq)
    chunk_iterator(missing_tms, upsert_chunk, label="Updating database...")

    tmp = Primer.select(Primer.seq).where((Primer.tm <= max_tm) &
                                          (Primer.tm >= min_tm))
    swga.message("{} primers have a melting temp within given range"
                 .format(tmp.count())) 

    # We mark all primers that pass criteria, then order by bg binding, then
    # select only the top few for location-finding
    valid_primers = (Primer
                     .update(active=True)
                     .where((Primer.seq << fgp) &
                            (Primer.seq << bgp) &
                            (Primer.seq << tmp))
                     .execute())

    swga.message("\nIdentified {} qualifying primers".format(valid_primers))
    
    # Order by bg binding, select up to num_primers, and reorder those by ratio
    subquery = (Primer
                .select()
                .where(Primer.active==True)
                .order_by(Primer.bg_freq, Primer.fg_freq.desc())
                .limit(num_primers))
    query = (Primer
             .select(Primer.seq)
             .where(Primer.seq << subquery)  # << means "IN"
             .order_by(Primer.ratio.desc()))

    # Un-mark the rest of the primers 
    Primer.update(active=False).where(~(Primer.seq << query)).execute()
    n_active = query.count()

    # Find binding sites in fg if they don't already have them
    if n_active > 0:
        missing_locs = list(Primer
                            .select()
                            .where((Primer.active == True) &
                                   (Primer.locations >> None))
                            .execute())
        if 0 < len(missing_locs) < 15:
            for p in progress.bar(missing_locs):
                p.locations = locate.primer_bind_sites(p.seq, fg_genome_fp)
        elif 0 < len(missing_locs):
            missing_locs = locate.primers_in_parallel(missing_locs,
                                                      fg_genome_fp)

        chunk_iterator(missing_locs, upsert_chunk,
                       label = "Updating primer binding sites...")

    swga.message("Marked {} primers as active (max specified as {})"
                 .format(n_active, num_primers))



