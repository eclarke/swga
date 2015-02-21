import json
import swga.database
from swga.database import Primer, update_in_chunks
from swga.melting import Tm
from swga.clint.textui import progress
import swga.locate as locate
import swga.primers
from swga.commands import Command

def main(argv, cfg_file):
    cmd = Command('activate', cfg_file=cfg_file)
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)

    primers = swga.primers.read_primer_list(
        cmd.input,
        cmd.fg_genome_fp,
        cmd.bg_genome_fp)

    activate_primers(primers, cmd.fg_genome_fp)


def activate_primers(primers, fg_genome_fp):
    update_tms(primers)
    update_locations(primers, fg_genome_fp)
    n_active = Primer.update(active=True).where(Primer.seq << primers).execute()
    swga.message("Marked {} primers as active.".format(n_active))
        

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
            p.locations = json.dumps(locate.binding_sites(p.seq, fg_genome_fp))
    elif 0 < len(primers):
        primers = locate.primers_in_parallel(primers,
                                             fg_genome_fp)
    update_in_chunks(primers, label = "Updating database... ")
