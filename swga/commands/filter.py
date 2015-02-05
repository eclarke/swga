import os
import json
import swga
from swga.commands import Command
from swga.primers import Primer, get_primer_locations
from swga.clint.textui import progress

def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_kmers(**cmd.args)


def filter_kmers(primer_db,
                 fg_genome_fp,
                 fg_length,
                 bg_length,
                 min_avg_fg_bind,
                 max_avg_bg_bind,
                 max_bg_bind,
                 min_tm,
                 max_tm,
                 num_primers):

    swga.primers.init_db(primer_db)

    fg_min_freq = float(min_avg_fg_bind) * fg_length
    bg_max_freq = float(max_avg_bg_bind) * bg_length

    assert fg_min_freq > 0
    
    Primer.update(active = False).execute()
    
    fgp = Primer.select(Primer.seq).where(Primer.fg_freq >= fg_min_freq)
    swga.message("{} primers bind foreground genome with avg rate >= {} sites/bp"
                 .format(fgp.count(), min_avg_fg_bind)) 

    bgp = Primer.select(Primer.seq).where(Primer.bg_freq <= bg_max_freq)
    swga.message("{} primers bind background genome with avg rate <= {} sites/bp"
                 .format(bgp.count(), max_avg_bg_bind))
    
    tmp = Primer.select(Primer.seq).where((Primer.tm <= max_tm) &
                                          (Primer.tm >= min_tm))
    swga.message("{} primers have a melting temp within given range"
                 .format(tmp.count())) 

    valid_primers = (Primer
                     .update(active=True)
                     .where((Primer.seq << fgp) &
                            (Primer.seq << bgp) &
                            (Primer.seq << tmp))
                     .execute())

    swga.message("\nIdentified {} qualifying primers".format(valid_primers))

    subquery = (Primer
                .select(Primer.seq)
                .where(Primer.active==True)
                .order_by(Primer.bg_freq)
                .limit(num_primers))

    query = (Primer
             .select()
             .where(Primer.seq << subquery)  # << means "IN"
             .order_by(Primer.ratio.desc()))
    
    selected_primers = [p.seq for p in query.execute()]

    # Reset "active" marks
    Primer.update(active = False).execute()

    n_active = (Primer
                .update(active = True)
                .where(Primer.seq << selected_primers).execute())
    if n_active > 0:
        empty_locations = Primer.select().where((Primer.active == True) &
                                                (Primer.locations >> None))

        for primer in progress.bar(empty_locations, 
                                   expected_size=empty_locations.count(),
                                   label="Finding fg binding sites: "):
            primer.locations = json.dumps(get_primer_locations(primer.seq,
                                                               fg_genome_fp))
            primer.save()

    swga.message("Marked {} primers as active (max specified as {})"
                 .format(n_active, num_primers))


