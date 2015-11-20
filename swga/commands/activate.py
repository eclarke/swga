from swga import warn
from swga.primers import Primers
from swga.commands._command import Command


class Activate(Command):

    def run(self):
        primers = Primers(self.input)

        try:
            (primers
                .update_melt_temps()
                .update_locations(self.fg_genome_fp)
                .activate())
        except AttributeError as e:
            warn("Error updating database: '{}'".format(e.message))
            raise e
