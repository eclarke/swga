import swga.database
import swga.primers
from swga import warn
from swga.primers import Primers
from swga.commands import Command


def main(argv, cfg_file):
    cmd = Command('activate', cfg_file)
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)

    primers = Primers(cmd.input)

    try:
        (primers
            .update_melt_temps()
            .update_locations(cmd.fg_genome_fp)
            .activate())
    except AttributeError as e:
        warn("Error updating database: '{}'".format(e.message))
        warn(
            "Sometimes this happens if you're using a database created with "
            "an older version of swga. Please try creating a new workspace "
            "with `swga init` in another folder and adding these primers to "
            "that database instead."
        )
        raise e
