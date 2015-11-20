from ConfigParser import SafeConfigParser
import abc
import click
import argutils
from swga import (
    __version__,
    warn,
    message
)
import swga.utils as utils


class Command(object):

    def __init__(self, name, cfg_file, metadata):
        """Create a new command with the given name and argument string.

        :param name: the command name, used to find the appropriate section in
        the config file
        :param cfg_file: the path to the config file
        :param metadata: a swga.metadata namedtuple generated from
        database.get_metadata()
        """
        spec = utils.specfile(name)
        opts = argutils.read.from_yaml(spec)

        self.name = name
        self.parser = argutils.export.to_argparser(name, opts)
        config = SafeConfigParser()
        if config.read(cfg_file) and config.has_section(name):
            self.parser = argutils.set_parser_defaults(self.parser, config)
        self._set_metadata(metadata)

    @abc.abstractmethod
    def run(self):
        return

    def _set_metadata(self, meta):
        self.primer_db = meta.db_name
        self.fg_genome_fp = meta.fg_file
        self.bg_genome_fp = meta.bg_file
        self.exclude_fp = meta.ex_file
        self.fg_length = meta.fg_length
        self.bg_length = meta.bg_length

    def parse_args(self, argv, quiet=False):
        args, unknown = self.parser.parse_known_args(argv)
        self.args = vars(args)
        self.unknown = unknown
        for k, v in self.args.items():
            setattr(self, k, v)
        if not quiet:
            self._pprint_args()
        if len(self.args) > 0 and all(v is None for k, v in self.args.items()):
            warn(
                "[swga %s]: All parameters are missing- ",
                "this may indicate a corrupt or missing parameters file."
                % self.name)

    def _pprint_args(self):
        message(
            click.style("swga {}, v{}".format(self.name, __version__),
                        fg='green'))
        for arg, val in self.args.iteritems():
            if val is None or val == "":
                val = click.style("None", fg='red')
            message(click.style("  - {}: {}".format(arg, val), fg='blue'))
