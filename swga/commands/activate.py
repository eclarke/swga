import swga.database
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

    swga.database.Primer.update_Tms(primers)
    swga.database.Primer.update_locations(primers, cmd.fg_genome_fp)
    n_activated = swga.database.Primer.activate(primers)
    swga.message("Marked {} primers as active.".format(n_activated))
