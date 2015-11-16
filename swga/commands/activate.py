from swga import warn
from swga.primers import Primers
from swga.commands._command import Command


class Activate(Command):

    def __init__(self, argv):
        super(Activate, self).__init__('activate')
        self.parse_args(argv)

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
